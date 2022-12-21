#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
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

import os.path

from nomad.metainfo import MSection, Package, Quantity, SubSection

from nomad.datamodel.data import ArchiveSection


m_package = Package(name='downloads')


class Download(MSection):
    url = Quantity(type=str, description='''
        A valid and downloadable URL. Will be downloaded on the servers that
        run this entries processing (e.g. NOMAD servers). The files will be
        added to the given output directory.''')

    output = Quantity(type=str, default='./', description='''
        A relative path that denotes the file to download the given URL to.
        Any parent directories will be created if they do not exist.
        Files that are marked to be extracted will be downloaded and extracted into
        the parent directory of the given file path.''')

    extract = Quantity(type=bool, description='''
        If the given URL denotes a compressed file and this flag is set to true,
        the downloaded file will be extracted and removed. Supported file formats
        are `.zip`, `.tgz`, `.tar.gz`.''')


class Downloads(ArchiveSection):
    '''
    Allows you to upload a very small file that will add very large files to your upload.
    Imagine there are large file resources in the internet (e.g. on a data sharing service)
    that you need to add to your upload. This way you do not need to download those large
    files first, just to upload them to NOMAD.

    When this section is processed, it will download files from given URLs, add
    them to the upload, and trigger processing for given mainfiles.
    '''

    description = Quantity(
        type=str, description='''Provides some additional description for these downloads.''')

    mainfiles = Quantity(type=str, shape=['*'], description='''
        A list of relative paths that denote mainfiles. These files are subjected
        to NOMAD processing after all files have been downloaded and potentially
        extracted.''')

    skip_download = Quantity(type=bool, description='''
        If true, the downloads will not be performed and no processing is triggered.
        If false, this will be changed to true by the processing after performing
        the downloads.''')

    downloads = SubSection(section=Download, repeats=True, description='''
        Defines URLs and how to download them.''')

    def normalize(self, archive, logger):
        super(Downloads, self).normalize(archive, logger)

        from nomad.datamodel import EntryArchive, ServerContext

        # Configure the entry type and name if this is the data section
        if archive.data == self:
            archive.metadata.entry_type = 'Downloads'
            archive.metadata.entry_name = archive.metadata.mainfile

        # Test if the download was executed already
        if isinstance(archive.m_context, ServerContext):
            # Read the old existing archive, find the same downloads section, and compare
            # for `skip_download`
            try:
                archive_reader = archive.m_context.upload_files.read_archive(archive.metadata.entry_id)
                old_archive = EntryArchive.m_from_dict(
                    archive_reader[archive.metadata.entry_id].to_dict())

                path = []
                current = self
                while current.m_parent:
                    path.append((current.m_parent_sub_section, current.m_parent_index))
                    current = current.m_parent

                current = old_archive
                for sub_section, index in reversed(path):
                    current = current.m_get_sub_section(sub_section, index)

                if current.skip_download:
                    logger.info('skipping downloads')
                    return
            except KeyError:
                # The old archive does not yet exist
                pass

        import pathlib
        import urllib.request
        import zipfile
        import tarfile
        from nomad.files import auto_decompress

        # download and extract files
        skip_download = True
        raw_path = archive.m_context.raw_path()
        for download in self.downloads:
            output = os.path.join(raw_path, download.output)
            download_logger = logger.bind(data=dict(
                path=output, abs_path=os.path.abspath(output)))

            if not os.path.abspath(output).startswith(os.path.abspath(raw_path)):
                download_logger.error('output path is not within the upload')
                continue

            if os.path.exists(output):
                download_logger.info('downloaded file does already exist')
                continue

            directory = os.path.dirname(output)
            if not os.path.exists(directory):
                pathlib.Path(directory).mkdir(parents=True)

            if not os.path.isdir(directory):
                download_logger.error('output path is not in a directory')
                continue

            try:
                urllib.request.urlretrieve(download.url, output)
            except Exception as e:
                download_logger.error('could not download url', exc_info=e)
                skip_download = False  # try again on next processing
                continue

            download_logger.info('downloaded file')

            if download.extract:
                format = auto_decompress(output)
                if format == 'zip':
                    with zipfile.ZipFile(output) as zf:
                        zf.extractall(directory)
                    download_logger.info('zip extracted file')
                elif format == 'tar':
                    with tarfile.open(output) as tf:
                        tf.extractall(directory)
                    download_logger.info('tar extracted file')
                else:
                    download_logger.error('unknown compression format')

                os.remove(output)

        # trigger processing for mainfiles
        mainfiles = self.mainfiles if self.mainfiles else []
        for mainfile in mainfiles:
            try:
                archive.m_context.process_updated_raw_file(mainfile, allow_modify=True)
            except Exception as e:
                logger.error('could not trigger processing', mainfile=mainfile, exc_info=e)
                skip_download = False
            else:
                logger.info('triggered processing', mainfile=mainfile)

        # skip download on next processing
        self.skip_download = skip_download


m_package.__init_metainfo__()
