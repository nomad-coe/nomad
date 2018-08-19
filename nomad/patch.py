import systax.analysis.symmetryanalyzer


def segfault_protect_patch(f, *args, **kwargs):
    return f(*args, **kwargs)

systax.analysis.symmetryanalyzer.segfault_protect = segfault_protect_patch
