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
from __future__ import annotations

import functools
import hashlib
import os.path
from multiprocessing import Pool, Lock, Manager
from typing import Iterable, Callable

from nomad.config import config
from nomad.archive import to_json, read_archive
from nomad.archive.storage_v2 import ArchiveWriter as ArchiveWriterNew
from nomad.files import StagingUploadFiles, PublicUploadFiles
from nomad.infrastructure import setup
from nomad.processing import Upload

lock = Lock()


def flush(*args, **kwargs):
    with lock:
        print(flush=True, *args, **kwargs)


class Counter:
    def __init__(self, total: int):
        manager = Manager()
        self.counter = manager.Value('i', 0)
        self.lock = manager.Lock()
        self.total = total

    def increment(self):
        with self.lock:
            self.counter.value += 1
            return f'[{self.counter.value}/{self.total}]'


def convert_archive(
    original_path: str,
    *,
    transform: Callable = None,
    overwrite: bool = False,
    delete_old: bool = False,
    counter: Counter = None,
    force_repack: bool = False,
):
    """
    Convert an archive of the old format to the new format.

    This code defines a function convert_archive that takes in an original file path, an
    optional transformation function, and a flag to indicate whether to overwrite
    existing files. The function checks if the file exists and is in the correct
    format, then reads the file and performs a conversion if necessary. Finally, it
    handles errors and overwrites the original file with the converted version if
    necessary.

    `transform` is a function that takes in the original file path as an argument, and
    returns the transformed file path. If `transform` is not provided, the original
    file path is overwritten with the converted version.

    `overwrite` is a boolean flag that specifies whether to overwrite the target
    file if it already exists. The default value is `False`.

    Args:
        original_path (str): The path to the original archive file.
        transform (Callable, optional): A function to transform the file name. Defaults to None.
        overwrite (bool, optional): Whether to overwrite existing files. Defaults to False.
        delete_old (bool, optional): Whether to delete the old file after conversion. Defaults to False.
        counter (Counter, optional): A counter to track the progress of the conversion. Defaults to None.
        force_repack (bool, optional): Force repacking the archive that is already in the new format. Defaults to False.
    """
    prefix: str = counter.increment() if counter else ''

    if not os.path.exists(original_path):
        flush(f'{prefix} [ERROR] File not found: {original_path}')
        return

    if not original_path.endswith('.msg'):
        flush(f'{prefix} [ERROR] Not a msgpack file: {original_path}')
        return

    original_path = os.path.abspath(original_path)

    if not force_repack:
        with open(original_path, 'rb') as f:
            magic_bytes = f.read(ArchiveWriterNew.magic_len)

        if magic_bytes == ArchiveWriterNew.magic:
            flush(
                f'{prefix} [INFO] Skipping as already in the new format: {original_path}'
            )
            return

    def safe_remove(path: str):
        if not path:
            return
        try:
            os.remove(path)
        except OSError:
            pass

    try:
        tmp_path = ''
        with read_archive(original_path, use_blocked_toc=False) as reader:
            flush(f'{prefix} [INFO] Converting: {original_path}')
            tmp_path = (
                f'{original_path}.{hashlib.md5(original_path.encode()).hexdigest()}'
            )

            with ArchiveWriterNew(
                tmp_path, len(reader), config.archive.toc_depth
            ) as writer:
                for uuid, entry in reader.items():
                    writer.add(uuid, to_json(entry))
    except Exception as e:
        flush(f'{prefix} [ERROR] Failed to convert {original_path}: {e}')
        safe_remove(tmp_path)
    else:
        new_path = transform(original_path) if transform else original_path
        if os.path.exists(new_path):
            if not overwrite:
                flush(f'{prefix} [ERROR] File already exists: {new_path}')
                safe_remove(tmp_path)
                return

            safe_remove(new_path)

        # the old path and the new path could be the same
        if delete_old and os.path.exists(original_path):
            safe_remove(original_path)

        os.rename(tmp_path, new_path)


def convert_folder(
    folders: str | list[str],
    *,
    processes: int = os.cpu_count(),
    transform: Callable = None,
    if_include: Callable = None,
    overwrite: bool = False,
    delete_old: bool = False,
    force_repack: bool = False,
):
    """
    Convert archives in the specified folder to the new format using parallel processing.

    transform(file_path:str) -> str

    if_include(file_path:str) -> bool

    Args:
        folders (str | list[str]): The path to the folder(s) containing the archives.
        processes (int): The number of parallel processes to use (default is 1).
        transform (Callable): A function to transform the file name (default is None).
        if_include (Callable): A function to filter the files to be converted (default is None).
        overwrite (bool): Whether to overwrite existing files (default is False).
        delete_old (bool): Whether to delete the old file after conversion (default is False).
        force_repack (bool): Force repacking the archive (default is False).
    """
    file_list: list = []

    if isinstance(folders, str):
        folders = [folders]

    if if_include is None:

        def _if_include(_):
            return True

        if_include = _if_include

    for folder in folders:
        if not os.path.exists(folder):
            flush(f'[ERROR] Folder not found: {folder}')
            continue

        for root, _, files in os.walk(folder):
            for file in files:
                full_name = os.path.join(root, file)
                if file.endswith('.msg') and if_include(full_name):
                    file_list.append(full_name)

    if not file_list:
        return

    counter = Counter(len(file_list))

    with Pool(processes=processes) as pool:
        pool.map(
            functools.partial(
                convert_archive,
                transform=transform,
                overwrite=overwrite,
                delete_old=delete_old,
                counter=counter,
                force_repack=force_repack,
            ),
            file_list,
        )


def convert_upload(
    uploads: Upload | Iterable[Upload] | str | Iterable[str],
    *,
    processes: int = os.cpu_count(),
    transform: Callable = None,
    if_include: Callable = None,
    overwrite: bool = False,
    delete_old: bool = False,
    force_repack: bool = False,
):
    """
    Function to convert an upload with the given upload_id to the new format.

    transform(file_path:str) -> str

    if_include(file_path:str) -> bool

    Args:
        uploads (str): The upload(s) to be converted.
        processes (int, optional): The number of processes to use for conversion. Defaults to 1.
        transform (Callable, optional): A function to apply transformation to the file name. Defaults to None.
        if_include (Callable, optional): A function to filter the files to be converted. Defaults to None.
        overwrite (bool, optional): Whether to overwrite existing files. Defaults to False.
        delete_old (bool, optional): Whether to delete the old file after conversion. Defaults to False.
        force_repack (bool, optional): Force repacking the existing archive (in new format). Defaults to False.
    """
    if isinstance(uploads, (str, Upload)):
        uploads = [uploads]

    if not uploads:
        return

    if isinstance(uploads[0], str):  # type: ignore
        setup()
        uploads = Upload.objects(upload_id__in=uploads)

    all_folders: list = []

    for upload in uploads:
        assert isinstance(upload, Upload)
        upload_class = PublicUploadFiles if upload.published else StagingUploadFiles
        base_folder = upload_class.base_folder_for(upload.upload_id)
        if not os.path.exists(base_folder):
            flush(f'[ERROR] Base folder not found for upload: {upload.upload_id}')
            continue

        all_folders.append(os.path.abspath(base_folder))

    convert_folder(
        all_folders,
        processes=processes,
        transform=transform,
        if_include=if_include,
        overwrite=overwrite,
        delete_old=delete_old,
        force_repack=force_repack,
    )


if __name__ == '__main__':

    def rename(path):
        return path.replace('v1.msg', 'v3.msg')

    def only_new(path):
        return 'old' not in path

    # filter files to be converted
    # convert_folder(
    #     '/home/theodore/PycharmProjects/nomad-FAIR/.volumes/fs/staging',
    #     transform=rename,
    #     if_include=only_new,
    # )
    #
    # # delete old archives after successful conversion
    # convert_folder(
    #     '/home/theodore/PycharmProjects/nomad-FAIR/.volumes/fs/staging',
    #     transform=rename,
    #     delete_old=True,
    # )
    #
    # # overwrite any existing files
    # convert_folder(
    #     '/home/theodore/PycharmProjects/nomad-FAIR/.volumes/fs/staging',
    #     transform=rename,
    #     overwrite=True,
    # )

    # overwrite existing archive
    convert_folder(
        '/home/theodore/PycharmProjects/nomad-FAIR/.volumes/fs/staging',
        overwrite=True,
    )
