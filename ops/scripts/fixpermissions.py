import multiprocessing as mp
import sys
import os
import os.path
import stat


if __name__ == '__main__':
    def fix_permissions(path):
        has_problems = False
        for item in os.listdir(path):
            filepath = os.path.join(path, item)
            if os.path.islink(filepath):
                continue
            stats = os.stat(filepath)
            has_problems |= not bool(stats.st_mode & stat.S_IROTH)
            if stats.st_mode & stat.S_IFDIR:
                has_problems |= not bool(stats.st_mode & stat.S_IXOTH)
        
        if has_problems:
            print('fixing problems for %s' % path)
            os.system('find %s -type d -exec chmod +x {} \\;' % path)
            os.system('find %s -exec chmod +r {} \\;' % path)

    if len(sys.argv) == 2:
        path = sys.argv[1]
        paths = [os.path.join(path, item) for item in os.listdir(path)]
    elif len(sys.argv) == 1:
        print('no path given')
    else:
        paths = sys.argv[1:]

    with mp.Pool(5) as p:
        p.map(fix_permissions, paths)
