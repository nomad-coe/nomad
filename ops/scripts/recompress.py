import os
import os.path
import click
import zipfile
import uuid
import shutil
import logging


@click.command()
@click.option('--src', default='/fairdi/migration/fs/migration_packages', type=str)
@click.option('--tmp', default='/scratch/fairdi/tmp/recompress', type=str)
def command(src, tmp):
    for root, _, files in os.walk(src):
        for file in files:
            try:
                path = os.path.join(root, file)
                if not path.endswith('.zip'):
                    continue

                tmp_path = os.path.join(tmp, str(uuid.uuid4()))
                new_path = tmp_path + '_new.zip'

                print(path)
                with zipfile.ZipFile(path, mode='r', allowZip64=True) as zip_file:
                    os.makedirs(tmp_path)
                    zip_file.extractall(tmp_path)

                with zipfile.ZipFile(
                        new_path, mode='w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
                    for tmp_root, _, tmp_files in os.walk(tmp_path):
                        for tmp_file in tmp_files:
                            zip_file.write(os.path.join(tmp_root, tmp_file), os.path.relpath(os.path.join(tmp_root, tmp_file), tmp_path))

                shutil.rmtree(tmp_path)
                shutil.move(new_path, path)
            except Exception as e:
                logging.error('exception', exc_info=e)
            finally:
                try:
                    shutil.rmtree(tmp_path)
                    os.remove(new_path)
                except Exception:
                    pass

            logging.info('recompressed ' + path)


if __name__ == '__main__':
    command()
