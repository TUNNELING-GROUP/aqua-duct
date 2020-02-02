# -*- coding: utf-8 -*-


from aquaduct.apps.data import GCS


################################################################################
# import or create memory decorator

if GCS.cachedir:
    from joblib import Memory
    memory_cache = Memory(cachedir=GCS.cachedir,
                          verbose=0)
    # mmap have to be switched off, otherwise smoothing does not work properly
    # memory_cache = Memory(cachedir=GCS.cachedir, mmap_mode='r', verbose=0)
    memory = memory_cache.cache
elif GCS.cachemem:
    pass
else:
    pass




################################################################################
# flag sandwich as imported
GCS.sandwich_import = True
################################################################################
