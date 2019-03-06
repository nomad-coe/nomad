import sys
from multiprocessing import Pool
import shutil

if __name__ == '__main__':
    paths = sys.argv[1:]
    with Pool(23) as p:
        p.map(shutil.rmtree, paths)
