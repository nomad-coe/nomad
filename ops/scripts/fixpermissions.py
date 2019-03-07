import multiprocessing as mp
import sys
import os
import os.path
import stat


if __name__ == '__main__':
    paths = sys.argv[1:]

    def fix_permissions(path):
        has_problems = False
        for item in os.listdir(path):
            filepath = os.path.join(path, item)
            stats = os.stat(filepath)
            has_problems |= not bool(stats.st_mode & stat.S_IROTH)
            if stats.st_mode & stat.S_IFDIR:
                has_problems |= not bool(stats.st_mode & stat.S_IXOTH)

        print('fixing problems for %s' % path)
        os.system('find %s -type d -exec chmod +x {} \\;' % path)
        os.system('find %s -exec chmod +r {} \\;' % path)

    with mp.Pool(5) as p:
        p.map(fix_permissions, paths)
