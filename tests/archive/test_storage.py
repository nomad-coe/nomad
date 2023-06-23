import gc
import os
import random
import string
import time
from io import BytesIO

import matplotlib.pyplot as plt
import pytest

from nomad.archive import write_archive
from nomad.archive.storage_v2 import ArchiveReadCounter, ArchiveReader, to_json, ArchiveList, ArchiveDict

# set matplotlib font size
plt.rcParams.update({'font.size': 6})


def generate_random_json(depth=10, width=4, simple=False):
    seed = random.random()

    if depth == 0 or (simple and seed < 0.4):
        return random.choice([
            random.randint(1, 100),
            random.random(),
            random.choice([True, False]),
            ''.join(random.choices(string.ascii_letters + string.digits, k=random.randint(5, 10)))])

    if seed < 0.7:
        obj = {}
        for _ in range(width):
            key = ''.join(random.choices(string.ascii_lowercase, k=random.randint(5, 10)))
            value = generate_random_json(depth - 1, width, True)
            obj[key] = value
        return obj

    lst = []
    for _ in range(width):
        value = generate_random_json(depth - 1, width, True)
        lst.append(value)
    return lst


def find_all_paths(json_obj, path=None, paths_list=None):
    if paths_list is None:
        paths_list = []

    if isinstance(json_obj, (dict, ArchiveDict)):
        for key, value in json_obj.items():
            new_path = [key] if not path else path + [key]
            find_all_paths(value, new_path, paths_list)
    elif isinstance(json_obj, (list, ArchiveList)):
        for index, value in enumerate(json_obj):
            new_path = [index] if not path else path + [index]
            find_all_paths(value, new_path, paths_list)

    if path:
        paths_list.append(path)

    return paths_list


def to_path(container, path):
    for key in path:
        container = container[key]
    return container


def test_generate_random_json(monkeypatch, tmp):
    monkeypatch.setattr('nomad.config.archive.use_new_writer', True)
    monkeypatch.setattr('nomad.config.archive.read_buffer_size', 4096)
    monkeypatch.setattr('nomad.config.archive.small_obj_optimization_threshold', 1024)

    depth = 6
    width = 6
    archive = {f'id': generate_random_json(depth, width)}
    path = find_all_paths(archive)

    def benchmark(_toc_depth):
        random.shuffle(path)
        f = BytesIO()
        write_archive(f, 1, [(f'id', archive)], _toc_depth)

        total_size = len(f.getvalue())

        counter = ArchiveReadCounter()

        with ArchiveReader(f, use_blocked_toc=False, counter=counter) as reader:
            accu = []
            for i in path:
                _ = to_path(reader['id'], i)
                assert to_path(reader['id'], i) == to_path(archive, i)
                accu.append(counter() / total_size)

            # assert counter() <= total_size

        counter.clear()

        with ArchiveReader(f, use_blocked_toc=False, counter=counter) as reader:
            _ = to_json(reader['id'])  # pylint: disable=no-member

        # assert counter() <= total_size

        return accu, total_size

    repeat = 10
    fig = plt.figure(figsize=(6, 3 * depth))
    for d in range(2, depth + 1):
        ax = fig.add_subplot(depth - 1, 1, d - 1)
        total_size = 0
        for _ in range(repeat):
            record, size = benchmark(d)
            x = [float(i) / len(record) for i in range(len(record))]
            ax.plot(x, record)
            total_size += size

        ax.set_xlabel('fraction of paths visited')
        ax.set_ylabel('fraction of bytes read')
        ax.set_title(f'archive size {total_size / repeat / 1024 / 1024:.2} MB d/w={depth}/{width}, toc depth={d - 1}')
        ax.set_xlim(-.01, .1)
        ax.set_ylim(0, 1.05)
        ax.grid(True)
    fig.tight_layout()
    plt.savefig(f'toc.png')


parent_folder = '/nomad/fairdi/tmp/thchang/archive_test/'
depth = 8
threshold = (4, 8, 16, 32, 64, 128, 256, 512)


@pytest.mark.skip
def test_folder_access(monkeypatch):
    archive = {f'id': generate_random_json(depth, 12)}
    paths = find_all_paths(archive)

    with open(f'{parent_folder}paths.txt', 'w') as f:
        for path in paths:
            f.write(','.join([str(v) for v in path]) + ',\n')

    from nomad.archive.storage_v2 import write_archive as write_archive_v2
    from nomad import config
    for j in threshold:
        monkeypatch.setattr('nomad.config.archive.small_obj_optimization_threshold', j * 1024)
        for i in range(2, depth):
            file_name = f'archive-{i}-{config.archive.small_obj_optimization_threshold}.msg'
            file_path = f'{parent_folder}{file_name}'
            write_archive_v2(file_path, [(f'id', archive)], i)


def measure(small_obj, toc_depth, paths):
    accu_time = []
    accu_size = []
    avg_read = []
    total_time = 0
    counter = ArchiveReadCounter()
    file_path = f'{parent_folder}archive-{toc_depth}-{small_obj}.msg'
    file_size = os.path.getsize(file_path)
    with ArchiveReader(file_path, use_blocked_toc=False, counter=counter) as reader:
        entry = reader['id']
        for path in paths:
            start = time.monotonic_ns()
            _ = to_json(to_path(entry, path))
            total_time += time.monotonic_ns() - start
            accu_time.append(total_time / 1e9)
            accu_size.append(counter() / file_size * 100)
            avg_read.append(counter.bytes_per_call() / 1024)
    gc.collect()
    return accu_time, accu_size, avg_read


@pytest.mark.skip
def test_read_archive(monkeypatch):
    monkeypatch.setattr('nomad.config.archive.fast_loading', True)

    paths = []
    with open(f'{parent_folder}paths.txt', 'r') as f:
        while line := f.readline():
            paths.append([(int(v) if v.isdigit() else v) for v in line.split(',') if v and v != '\n'])

    rows = len(threshold)
    cols = depth - 2

    fig = plt.figure(figsize=(4 * cols, 3 * rows), dpi=300)

    pic_pool = {}
    counter = 0
    for j in range(rows):
        small_obj = threshold[j] * 1024
        for i in range(depth - 2):
            toc_depth = i + 2
            counter += 1
            ax1 = fig.add_subplot(rows, cols, counter)
            ax2 = ax1.twinx()
            ax3 = ax1.twinx()
            ax1.set_xlabel('number of paths visited')
            ax1.set_ylabel('accu. time (s)')
            ax2.set_ylabel('accu. size (%)')
            ax3.set_ylabel('KB per read call')
            ax1.set_title(f'toc={toc_depth}, small obj={threshold[j]} k')
            ax1.grid(True)
            ax1.yaxis.label.set_color('blue')
            ax2.yaxis.label.set_color('red')
            ax3.yaxis.label.set_color('green')
            ax3.spines.right.set_position(("axes", 1.2))
            pic_pool[f'archive-{toc_depth}-{small_obj}'] = (ax1, ax2, ax3)

    samples = 100_000
    x = [i + 1 for i in range(samples)]

    for _ in range(20):
        for i in range(depth - 2):
            toc_depth = i + 2
            for j in range(rows):
                small_obj = threshold[j] * 1024
                ax1, ax2, ax3 = pic_pool[f'archive-{toc_depth}-{small_obj}']
                accu_time, accu_size, avg_bytes = measure(small_obj, toc_depth, random.choices(paths, k=samples))
                ax1.semilogx(x, accu_time, color='blue')
                ax2.semilogx(x, accu_size, color='red')
                ax3.semilogx(x, avg_bytes, color='green')

        fig.tight_layout()
        plt.savefig(f'{parent_folder}toc_archive.png')
