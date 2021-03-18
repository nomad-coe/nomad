#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

'''
Functionality for "lazy loading". A lazy loaded module is registered in sys.modules as
loaded, but is not actually run intil needed, i.e. when an attribute which requires the
module to be imported is accessed. This is useful for optimization purposes, but needs to
be used with caution.

Since the actual import of a lazy loaded module is postponend until it is actually needed,
exceptions occuring during the actual import can happen in a location where you may not
have anticipated it. Threfore, exceptions raised during the lazy import are wrapped in a
dedicated exception class, LazyImportError, to avoid potential import exceptions being
misclassified.
'''

import sys
from types import ModuleType
from typing import Set
from importlib._bootstrap import _ImportLockContext
import importlib
from nomad import config


_not_yet_imported_lazy_module_names: Set[str] = set()


class LazyImportError(Exception):
    pass


class _LazyModule(ModuleType):
    '''
    Base class for lazy modules.
    '''
    def __getattribute__(self, attr):
        '''
        Overrides the standard method, to trigger the actual import of the lazy module
        when a non-trivial attribute is accessed.

        Note: when we actually import the module, we will replace this method with the
        standard method, thereby reverting to standard behaviour.
        '''
        if attr not in ('__name__', '__class__', '__spec__', '__repr__', '__file__'):
            try:
                # In case the attribute we are trying to access is actually another
                # lazy-loaded module, just return it.
                name = '%s.%s' % (self.__name__, attr)
                return sys.modules[name]
            except KeyError:
                pass
            # No, we have to actually load the module now!
            _actually_import(self.__name__)
        return ModuleType.__getattribute__(self, attr)  # Standard functionality

    def __setattr__(self, attr, value):
        '''
        Overrides the standard method, to trigger the actual import of the lazy module
        when any attribute is set on the module.

        Note: when we actually import the module, we will replace this method with the
        standard method, thereby reverting to standard behaviour.
        '''
        _actually_import(self.__name__)
        return ModuleType.__setattr__(self, attr, value)  # Standard functionality


def _actually_import(module_name):
    '''
    Actually import the lazy module. Also make sure that all its parent modules are
    imported if needed - and they need to be imported in the right order.
    '''
    with _ImportLockContext():
        parts = module_name.split('.')
        base_name = ''
        for part in parts:
            if base_name:
                base_name += '.'
            base_name += part
            if base_name in _not_yet_imported_lazy_module_names:
                # This level is a lazy module, and it has not yet been loaded. Load it!
                _not_yet_imported_lazy_module_names.remove(base_name)
                module = sys.modules[base_name]
                # Restore __getattribute__ and __setattr__ to original functionality
                module_class = type(module)
                module_class.__getattribute__ = ModuleType.__getattribute__
                module_class.__setattr__ = ModuleType.__setattr__
                # Remove the fake __file__ attribute we set initially
                ModuleType.__delattr__(module, '__file__')
                # Actually import the module
                try:
                    importlib.reload(module)
                except Exception as e:
                    # Wrap the exception, to avoid potential exception misclassification.
                    err_msg = f'Error occured during loading of lazy module {base_name}: {e}'
                    raise LazyImportError(err_msg)


def _create_lazy_module(module_name):
    '''
    Create a dedicated class and instantiate it, and adds it to sys.modules
    '''
    class _LazyModuleSubclass(_LazyModule):
        def __repr__(self):
            return 'Lazily-loaded module %s' % self.__name__
    module = _LazyModuleSubclass(module_name)
    ModuleType.__setattr__(module, '__file__', None)
    sys.modules[module_name] = module
    _not_yet_imported_lazy_module_names.add(module_name)
    return module


def lazy_module(module_name):
    '''
    Call this to "lazily" import a module. Subsequent calls to import will succeed
    immediately, without the module actually being imported. The module is imported
    first when it is actually used (by accessing an attribute which requires the
    module to really be imported).

    The lazy import functionality can also be disabled using the setting
        nomad.config.enable_lazy_import = False
    When disabled, this method does nothing, and modules are imported "as usual".
    '''
    if not config.enable_lazy_import:
        return

    if module_name not in sys.modules:
        # Create a lazy module object and add it, without really loading it.
        # Also add a lazy module object for all parent modules, if needed.
        with _ImportLockContext():
            module = _create_lazy_module(module_name)
            while True:
                parent_module_name, _, submodule_name = module_name.rpartition('.')
                if not parent_module_name:
                    break
                # Fetch or create parent_module
                if parent_module_name in sys.modules:
                    parent_module = sys.modules[parent_module_name]
                    parent_was_already_created = True
                else:
                    parent_module = _create_lazy_module(parent_module_name)
                    parent_was_already_created = False
                # Set module as an attribute on the parent_module
                ModuleType.__setattr__(parent_module, submodule_name, module)
                if parent_was_already_created:
                    break
                # Parent had to be lazy loaded too -> we need to also check parent's parent.
                module = parent_module
                module_name = parent_module_name
    return
