# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from typing import cast, Any, Tuple, List, Set, Dict, Iterable
import os
import json
from datetime import datetime, timedelta
from packaging import version

from nomad import config, utils, datamodel, search
from nomad.config.models import BundleImportSettings, BundleExportSettings
from nomad.files import (
    zipfile, PathObject, UploadFiles, PublicUploadFiles, StagingUploadFiles,
    FileSource, BrowsableFileSource, CombinedFileSource, StreamedFileSource, DiskFileSource, ZipFileSource,
    json_to_streamed_file, bundle_info_filename, StandardJSONDecoder)
from nomad.processing.base import ProcessStatus
from nomad.processing.data import Upload, Entry, mongo_entry_metadata
from fastapi import HTTPException, status


class BundleExporter:
    def __init__(
            self, upload: Upload, export_as_stream: bool, export_path: str, zipped: bool, overwrite: bool,
            export_settings: BundleExportSettings):
        '''
        Class for exporting an upload as a *bundle*. Bundles are used to export and import
        uploads between different NOMAD installations. After instantiating a BundleExporter,
        use the `export_bundle` method to do the actual exporting.

        Arguments:
            upload:
                The upload to export
            export_as_stream: If the bundle should be exported as a stream, rather than saved
                to a file or folder. You must specify either `export_as_stream = True` or `export_path`.
                Further, `zipped` must be set to True when exporting as a stream.
                The stream is returned by the `export_bundle` function.
            export_path: Defines the output path, when exporting to disk rather than as a stream.
                You must specify either `export_as_stream = True` or `export_path`. Set to
                None if exporting as a stream.
            zipped: if the bundle should be zipped. Set to False to export the bundle to disk
                as an uncompressed folder. If exporting as a stream, zipped must be set to True.
            overwrite:
                If the target file/folder should be overwritten by this operation. Not
                applicable if `export_as_stream` is True.
            export_settings:
                Settings for controlling the bundle content. See the
                `config.BundeExportSettings` for applicable options.
                NOTE: the dictionary must specify a *complete* set of options.
        '''
        BundleExporter.check_export_settings(export_settings)
        self.upload = upload
        self.export_as_stream = export_as_stream
        self.export_path = export_path
        self.zipped = zipped
        self.overwrite = overwrite
        self.export_settings = export_settings

    @classmethod
    def check_export_settings(cls, export_settings: BundleExportSettings):
        assert export_settings.include_archive_files or export_settings.include_raw_files, (
            'Export must include the archive files or the raw files, or both')

    def export_bundle(self) -> Iterable[bytes]:
        # Safety checks
        if self.export_as_stream:
            assert self.export_path is None, 'Cannot have `export_path` set when exporting as a stream.'
            assert self.zipped, 'Must have `zipped` set to True when exporting as stream.'
        else:
            assert self.export_path is not None, 'You must specify either `export_as_stream = True` or `export_path`.'
            assert self.overwrite or not os.path.exists(self.export_path), '`export_path` alredy exists.'
        assert not self.upload.process_running or self.upload.current_process == 'publish_externally', (
            'Upload is being processed.')

        file_source = CombinedFileSource(self._get_file_sources())

        # Export
        if self.export_as_stream:
            return file_source.to_zipstream()
        else:
            # Create parent dir if it does not exist
            parent_dir = os.path.dirname(os.path.abspath(self.export_path))
            if not os.path.exists(parent_dir):
                os.makedirs(parent_dir)
            if self.zipped:
                file_source.to_zipfile(self.export_path, self.overwrite)
            else:
                file_source.to_disk(self.export_path, False, self.overwrite)
        return None

    def _get_file_sources(self) -> Iterable[FileSource]:
        ''' Generator which yields all the `FileSource` objects to pack in the bundle. '''
        # 1. The bundle info json
        bundle_info = self._create_bundle_info()
        yield StreamedFileSource(json_to_streamed_file(bundle_info, bundle_info_filename))

        # 2. Files from the upload dir
        for file_source in self.upload.upload_files.files_to_bundle(self.export_settings):
            yield file_source

    def _create_bundle_info(self):
        ''' Create the bundle_info.json data '''
        bundle_info: Dict[str, Any] = dict(
            upload_id=self.upload.upload_id,
            source=config.meta.dict(),  # Information about the source system, i.e. this NOMAD installation
            export_settings=self.export_settings.dict(),
            upload=self.upload.to_mongo().to_dict(),
            entries=[entry.to_mongo().to_dict() for entry in self.upload.successful_entries])
        # Handle datasets
        dataset_ids: Set[str] = set()
        for entry_dict in bundle_info['entries']:
            entry_datasets = entry_dict.get('datasets')
            if entry_datasets:
                if not self.export_settings.include_datasets:
                    entry_dict['datasets'] = None
                else:
                    dataset_ids.update(entry_datasets)
        if self.export_settings.include_datasets:
            bundle_info['datasets'] = [
                datamodel.Dataset.m_def.a_mongo.get(dataset_id=dataset_id).m_to_dict()
                for dataset_id in sorted(dataset_ids)]

        return bundle_info


class BundleImporter:
    def __init__(self, user: datamodel.User, import_settings: BundleImportSettings, embargo_length: int = None):
        '''
        Class for importing an upload from a *bundle*.

        Arguments:
            user:
                The user requesting the import. Used to check permissions. Of omitted, no
                permission checks are done.
            import_settings:
                Settings for controlling the bundle content. See the
                `BundleImportSettings` for applicable options.
                NOTE: the dictionary must specify a complete set of options.
            embargo_length:
                Used to set the embargo length. If set to None, the value will be imported
                from the bundle. The value should be between 0 and 36. A value of 0 means
                no embargo.
        '''
        self.user = user
        self.import_settings = import_settings
        self.embargo_length = embargo_length
        # Internals
        self.bundle_path: str = None
        self.bundle: BrowsableFileSource = None
        self.upload: Upload = None
        self.upload_files: UploadFiles = None
        self._bundle_info: Dict[str, Any] = None

    @classmethod
    def looks_like_a_bundle(cls, path):
        ''' Fast method to make a (very shallow) check if the object specified by `path` looks like a valid bundle. '''
        assert os.path.exists(path), f'Path not found: {path}'
        if os.path.isfile(path):
            try:
                if not path.lower().endswith('.zip'):
                    return False
                zip_file = zipfile.ZipFile(path, 'r')
                zip_file.getinfo(bundle_info_filename)
                return True
            except Exception:
                return False
        return os.path.isfile(os.path.join(path, bundle_info_filename))

    def check_api_permissions(self):
        '''
        Checks if the specified user is allowed to import a bundle via the api. Raises a
        HTTPException if not. This is a quick check, which does not require the bundle to be opened.
        '''
        if not self.user:
            return  # No permission checks

        is_admin = self.user.is_admin
        is_oasis = not is_admin and self.user.is_oasis_admin and config.bundle_import.allow_bundles_from_oasis

        if not is_admin and not is_oasis:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED, detail='User not authorized to import bundles')

        if not is_admin:
            for k, v in self.import_settings.dict().items():
                if v != config.bundle_import.default_settings.dict().get(k):
                    raise HTTPException(
                        status_code=status.HTTP_401_UNAUTHORIZED,
                        detail=f'Changing the setting {k} requires an admin user')

    def open(self, bundle_path: str):
        self.bundle_path = bundle_path
        if os.path.isdir(bundle_path):
            self.bundle = DiskFileSource(bundle_path)
        else:
            assert zipfile.is_zipfile(bundle_path), '`path` must define a folder or a zipfile.'
            zip_file = zipfile.ZipFile(bundle_path, 'r')
            self.bundle = ZipFileSource(zip_file)

    def create_upload_skeleton(self) -> Upload:
        '''
        Creates an upload "skeleton" to be populated with the imported data.
        '''
        # Some basic sanity checks
        is_admin = not self.user or self.user.is_admin
        is_oasis = not is_admin and self.user.is_oasis_admin and config.bundle_import.allow_bundles_from_oasis

        if is_oasis and not config.bundle_import.allow_unpublished_bundles_from_oasis:
            if not self.bundle_info.get('upload', {}).get('publish_time'):
                raise HTTPException(
                    status_code=status.HTTP_400_BAD_REQUEST,
                    detail=f'Bundles uploaded from an oasis must be published in the oasis first.')

        keys_exist(self.bundle_info, ('upload_id', 'upload.main_author'), 'Missing key in bundle_info.json: {key}')
        upload_id = self.bundle_info['upload_id']
        main_author = self.bundle_info['upload']['main_author']
        try:
            Upload.get(upload_id)
            assert False, f'Upload with id {upload_id} already exists'
        except KeyError:
            pass
        main_author_user = datamodel.User.get(user_id=main_author)
        assert main_author_user is not None, f'Invalid main_author: {main_author}'
        # Create the Upload
        self.upload = Upload.create(upload_id=upload_id, main_author=main_author_user)
        return self.upload

    def import_bundle(self, upload: Upload, running_locally) -> str:
        '''
        Import the data to the provided upload skeleton (the most costly step of the
        import process). Should only be invoked from a @process or a @process_local method on
        Upload (to get proper handling of process_status, errors, blocking other processes, etc).
        '''
        self.upload = upload
        logger = self.upload.get_logger(bundle_path=self.bundle_path)
        current_time = datetime.utcnow()
        new_datasets: List[datamodel.Dataset] = []
        dataset_id_mapping: Dict[str, str] = {}
        entry_data_to_index: List[datamodel.EntryArchive] = []  # Data to index in ES
        try:
            self._check_bundle_and_settings(running_locally)
            self._import_upload_mongo_data(current_time)
            if self.import_settings.include_datasets:
                new_datasets, dataset_id_mapping = self._import_datasets()
            entries = self._import_entries_mongo_data(current_time, dataset_id_mapping)
            self._import_files()
            self.bundle.close()
            entry_data_to_index = self._get_entry_data_to_index(entries)
            # Everything looks good - save to mongo.
            self.upload.save()
            for entry in entries:
                entry.save()
            # Index ES
            self._index_search(entry_data_to_index)

            if self.import_settings.delete_bundle_on_success:
                self.delete_bundle()

            # Possibly reprocess
            if self.import_settings.trigger_processing:
                return self._reprocess_upload()
            return None
        except Exception as e:
            logger.error('could not import bundle', exc_info=e)
            self.bundle.close()
            if self.import_settings.delete_bundle_on_fail:
                self.delete_bundle()
            if self.import_settings.delete_upload_on_fail:
                # Delete everything
                self.upload.delete_upload_local()  # Will also delete files, entries and remove from elastic search
                if new_datasets:
                    for dataset in new_datasets:
                        dataset.a_mongo.delete()
                return ProcessStatus.DELETED  # Used by the framework when running as a @process
            else:
                # Just ensure the upload is deleted from search
                with utils.timer(logger, 'upload deleted from index'):
                    search.delete_upload(self.upload.upload_id, refresh=True)
                raise e

    def close(self):
        if self.bundle:
            self.bundle.close()

    def delete_bundle(self):
        ''' Deletes the bundle file, and optionally, it's parent folder (if empty). '''
        if os.path.exists(self.bundle_path):
            PathObject(self.bundle_path).delete()
        if self.import_settings.delete_bundle_include_parent_folder:
            parent_folder = os.path.dirname(self.bundle_path)
            if not os.listdir(parent_folder):
                PathObject(parent_folder).delete()

    @property
    def bundle_info(self):
        if not self._bundle_info:
            assert self.bundle, 'Must open a bundle before getting the bundle info.'
            with self.bundle.open(bundle_info_filename, 'rt') as f:
                self._bundle_info = json.load(f, cls=StandardJSONDecoder)
        return self._bundle_info

    def _check_bundle_and_settings(self, running_locally: bool):
        ''' Perform various initial sanity checks. '''
        # Sanity checks of the settings
        assert not (self.import_settings.trigger_processing and running_locally), (
            'Cannot use `trigger_processing` when running locally.')
        # Sanity checks of the bundle
        required_keys_root_level = (
            'upload_id', 'source.version', 'source.commit', 'source.deployment', 'source.deployment_url',
            'export_settings.include_raw_files',
            'export_settings.include_archive_files',
            'export_settings.include_datasets',
            'upload._id', 'upload.main_author',
            'upload.upload_create_time', 'upload.process_status', 'upload.license',
            'upload.embargo_length',
            'entries')

        keys_exist(self.bundle_info, required_keys_root_level, 'Missing key in bundle_info.json: {key}')

        # Check version
        try:
            bundle_nomad_version = self.bundle_info['source']['version']
            assert version.parse(bundle_nomad_version) >= version.parse(config.bundle_import.required_nomad_version), (
                f'Bundle created in NOMAD version {bundle_nomad_version}, '
                f'required at least {config.bundle_import.required_nomad_version}')
        except Exception:
            assert False, 'Bad bundle version'

    def _import_upload_mongo_data(self, current_time):
        upload_dict = self.bundle_info['upload']
        assert self.upload.upload_id == self.bundle_info['upload_id'] == upload_dict['_id'], (
            'Inconsisten upload id information')
        published = upload_dict.get('publish_time') is not None
        if published:
            assert self.bundle_info['entries'], 'Upload published but no entries in bundle_info.json'
        # Check user references
        check_user_ids([upload_dict['main_author']], 'Invalid main_author: {id}')
        check_user_ids(upload_dict.get('coauthors', []), 'Invalid coauthor reference: {id}')
        check_user_ids(upload_dict.get('reviewers', []), 'Invalid reviewers reference: {id}')
        # Define which keys we think okay to copy from the bundle
        upload_keys_to_copy = [
            'upload_name', 'main_author', 'coauthors', 'reviewers', 'embargo_length', 'license',
            'from_oasis', 'oasis_deployment_url']
        if self.import_settings.keep_original_timestamps:
            upload_keys_to_copy.extend(('upload_create_time', 'publish_time',))
        try:
            # Update the upload with data from the json, and validate it
            update = {k: upload_dict[k] for k in upload_keys_to_copy if k in upload_dict}
            self.upload.modify(**update)
            self.upload.validate()
        except Exception as e:
            assert False, 'Bad upload json data: ' + str(e)
        # Manage timestamps
        current_time_plus_tolerance = current_time + timedelta(minutes=2)
        if published and not self.import_settings.keep_original_timestamps:
            self.upload.publish_time = current_time
        for timestamp in (
                self.upload.upload_create_time,
                self.upload.last_update,
                self.upload.complete_time,
                self.upload.publish_time):
            assert timestamp is None or timestamp < current_time_plus_tolerance, (
                'Timestamp is in the future')
        # Manage source info
        if self.import_settings.set_from_oasis:
            self.upload.from_oasis = True
            source_deployment_url = self.bundle_info['source']['deployment_url']
            assert source_deployment_url, 'No source deployment_url defined'
            if not self.upload.oasis_deployment_url:
                # Note, if oasis_deployment_url is set in the bundle_info, we keep this
                # value as it is, since it indicates that the upload has been imported from
                # somewhere else originally (i.e. source_deployment_url would not be the
                # original source)
                self.upload.oasis_deployment_url = source_deployment_url
        # Validate embargo settings
        if self.embargo_length is not None:
            self.upload.embargo_length = self.embargo_length  # Importing with different embargo
        assert type(self.upload.embargo_length) == int and 0 <= self.upload.embargo_length <= 36, (
            'Invalid embargo_length, must be between 0 and 36 months')

    def _import_datasets(self) -> Tuple[List[datamodel.Dataset], Dict[str, str]]:
        ''' Creates datasets from the bundle. '''
        required_keys_datasets = (
            'dataset_id', 'dataset_name', 'user_id')

        assert 'datasets' in self.bundle_info, 'Missing datasets definition in bundle_info.json'
        datasets = self.bundle_info['datasets']
        new_datasets: List[datamodel.Dataset] = []
        dataset_id_mapping: Dict[str, str] = {}  # Map from old to new id (usually the same)
        for dataset_dict in datasets:
            keys_exist(dataset_dict, required_keys_datasets, 'Missing key in dataset definition: {key}')
            check_user_ids([dataset_dict['user_id']], 'Invalid dataset creator id: {id}')
            dataset_id = dataset_dict['dataset_id']
            try:
                existing_dataset = datamodel.Dataset.m_def.a_mongo.get(
                    user_id=dataset_dict['user_id'],
                    dataset_name=dataset_dict['dataset_name'])
                # Dataset by the given dataset_name and user_id already exists
                dataset_id_mapping[dataset_id] = existing_dataset.dataset_id
                # Note, it may be that a dataset with the same dataset_name and creator
                # is created in both environments. In that case, we consider them
                # to be the "same" dataset, even if they do not have the same dataset_id.
                # Thus, in that case the dataset id needs to be translated.
                assert not existing_dataset.doi, (
                    f'Matched dataset {existing_dataset.dataset_id} has a DOI, cannot be updated')
            except KeyError:
                # Completely new dataset, create it
                new_dataset = datamodel.Dataset(**dataset_dict)
                new_dataset.a_mongo.save()
                new_datasets.append(new_dataset)
                dataset_id_mapping[dataset_id] = dataset_id

        return new_datasets, dataset_id_mapping

    def _import_entries_mongo_data(self, current_time, dataset_id_mapping) -> List[Entry]:
        ''' Creates mongo entries from the data in the bundle_info '''
        required_keys_entry_level = (
            '_id', 'upload_id', 'mainfile', 'parser_name', 'process_status', 'entry_create_time')

        entries = []
        for entry_dict in self.bundle_info['entries']:
            keys_exist(entry_dict, required_keys_entry_level, 'Missing key for entry: {key}')
            assert entry_dict['process_status'] in ProcessStatus.STATUSES_NOT_PROCESSING, (
                'Invalid entry `process_status`')
            # Check referential consistency
            assert entry_dict['upload_id'] == self.upload.upload_id, (
                'Mismatching upload_id in entry definition')
            expected_entry_id = utils.generate_entry_id(
                self.upload.upload_id, entry_dict['mainfile'], entry_dict.get('mainfile_key'))
            assert entry_dict['_id'] == expected_entry_id, (
                'Provided entry id does not match generated value')
            check_user_ids(entry_dict.get('entry_coauthors', []), 'Invalid entry_coauthor reference: {id}')

            # Instantiate an entry object from the json, and validate it
            entry_keys_to_copy = list(mongo_entry_metadata)
            entry_keys_to_copy.extend((
                'upload_id', 'errors', 'warnings', 'last_status_message',
                'current_process', 'complete_time', 'worker_hostname', 'celery_task_id'))
            try:
                update = {k: entry_dict[k] for k in entry_keys_to_copy if k in entry_dict}
                update['entry_id'] = entry_dict['_id']
                if not self.import_settings.keep_original_timestamps:
                    update['entry_create_time'] = current_time
                entry: Entry = Entry.create(**update)
                entry.process_status = entry_dict['process_status']
                entry.validate()
            except Exception as e:
                assert False, 'Bad entry json data: ' + str(e)
            # Instantiate an EntryMetadata object to validate the format
            try:
                if self.import_settings.include_datasets:
                    entry_datasets = entry_dict.get('datasets')
                    if entry_datasets:
                        entry.datasets = [
                            dataset_id_mapping[id] for id in entry_datasets] or None
                else:
                    entry.datasets = None
                entry.mongo_metadata(self.upload)
            except Exception as e:
                assert False, 'Invalid entry metadata: ' + str(e)
            entries.append(entry)
        return entries

    def _import_files(self):
        try:
            cls = PublicUploadFiles if self.upload.published else StagingUploadFiles
            assert not os.path.exists(cls.base_folder_for(self.upload.upload_id)), 'Upload folder already exists'
            self.upload_files = cls(self.upload.upload_id, create=True)

            for file_source in self.upload_files.files_from_bundle(self.bundle, self.import_settings):
                file_source.to_disk(self.upload_files.os_path, overwrite=True)

            if self.upload.published and self.embargo_length is not None:
                # Repack the upload
                PublicUploadFiles(self.upload.upload_id).re_pack(with_embargo=self.embargo_length > 0)
        except Exception:
            # Something went wrong. Delete the files and re-raise the original exception
            if self.upload_files:
                self.upload_files.delete()
            raise

    def _get_entry_data_to_index(self, entries: List[Entry]) -> List[datamodel.EntryArchive]:
        entry_data_to_index = []
        if self.import_settings.include_archive_files:
            for entry in entries:
                try:
                    entry_metadata = entry.full_entry_metadata(self.upload)
                    entry_data_to_index.append(
                        cast(datamodel.EntryArchive, entry_metadata.m_parent))
                except Exception as e:
                    assert False, 'Invalid metadata in archive entry: ' + str(e)
            self.upload_files.close()  # Because full_entry_metadata reads the archive files.
        return entry_data_to_index

    def _index_search(self, entry_data_to_index: List[datamodel.EntryArchive]):
        # Index in elastic search
        if entry_data_to_index:
            search.index(
                entry_data_to_index, update_materials=config.process.index_materials,
                refresh=True)

    def _reprocess_upload(self):
        return self.upload._process_upload_local(reprocess_settings=self.import_settings.process_settings)


def keys_exist(data: Dict[str, Any], required_keys: Iterable[str], error_message: str):
    '''
    Checks if the specified keys exist in the provided dictionary structure `data`.
    Supports dot-notation to access subkeys.
    '''
    for key in required_keys:
        current = data
        for sub_key in key.split('.'):
            assert sub_key in current, error_message.replace('{key}', key)
            current = current[sub_key]


def check_user_ids(user_ids: Iterable[str], error_message: str):
    '''
    Checks if all user_ids provided in the Iterable `user_ids` are valid. If not, raises an
    AssertionError with the specified error message. The string {id} in `error_message` is
    replaced with the bad value.
    '''
    for user_id in user_ids:
        user = datamodel.User.get(user_id=user_id)
        assert user is not None, error_message.replace('{id}', user_id)
