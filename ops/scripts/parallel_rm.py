import os
import sys
from multiprocessing import Pool

if __name__ == '__main__':
    def rm(path):
        os.system('rm -rf %s' % path)

    paths = sys.argv[1:]
    with Pool(23) as p:
        p.map(rm, paths)
