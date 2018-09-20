import systax.analysis.symmetryanalyzer


# A patch for the segfault protection of systax (internally uses protection for spglib calls.)
# We basically disable the protection. The multiprocessing based original protection.
# somehow interfers with the celery work infrastructure and leads to a deaklock. Its a TODO.
def segfault_protect_patch(f, *args, **kwargs):
    return f(*args, **kwargs)


systax.analysis.symmetryanalyzer.segfault_protect = segfault_protect_patch
