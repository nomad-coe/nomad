from nomad.files import UploadFiles
from nomad.app.api.auth import create_authorization_predicate
from nomad.archive_library.filedb import ArchiveFileDB


def get_dbs(upload_id):
    upload_files = UploadFiles.get(upload_id, create_authorization_predicate(upload_id))

    if upload_files is None:
        return []

    files = upload_files.archive_file_msg('X')
    msgdbs = [ArchiveFileDB(f) for f in files if f is not None]
    return msgdbs
