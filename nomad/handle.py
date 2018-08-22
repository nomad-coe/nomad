import logging

from nomad.processing import handle_uploads

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    handle_uploads()
