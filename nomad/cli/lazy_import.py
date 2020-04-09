from types import ModuleType
import sys
from importlib._bootstrap import _ImportLockContext
from six import raise_from
from importlib import reload as reload_module


__all__ = ['lazy_module', 'LazyModule', '_MSG']


_CLS_ATTRS = (
    '_lazy_import_error_strings', '_lazy_import_error_msgs', '_lazy_import_callables',
    '_lazy_import_submodules', '__repr__'
)

_DICT_DELETION = ('_lazy_import_submodules',)

_MSG = ("{caller} attempted to use a functionality that requires module "
        "{module}, but it couldn't be loaded. Please install {install_name} "
        "and retry.")


class LazyModule(ModuleType):
    def __getattribute__(self, attr):
        if attr not in ('__name__', '__class__', '__spec__'):
            try:
                name = '%s.%s' % (self.__name__, attr)
                return sys.modules[name]
            except KeyError:
                pass

            try:
                return type(self)._lazy_import_callables[attr]
            except (AttributeError, KeyError):
                _load_module(self)
        return super(LazyModule, self).__getattribute__(attr)

    def __setattr__(self, attr, value):
        _load_module(self)
        return super(LazyModule, self).__setattr__(attr, value)


def _clean_lazy_submodule_refs(module):
    module_class = type(module)
    for entry in _DICT_DELETION:
        try:
            names = getattr(module_class, entry)
        except AttributeError:
            continue
        for name in names:
            try:
                super(LazyModule, module).__delattr__(name)
            except AttributeError:
                pass


def _clean_lazymodule(module):
    module_class = type(module)
    _clean_lazy_submodule_refs(module)

    module_class.__getattribute__ = ModuleType.__getattribute__
    module_class.__setattr__ = ModuleType.__setattr__
    class_attrs = {}
    for attr in _CLS_ATTRS:
        try:
            class_attrs[attr] = getattr(module_class, attr)
            delattr(module_class, attr)
        except AttributeError:
            pass
    return class_attrs


def _reset_lazy_submodule_refs(module):
    module_class = type(module)
    for entry in _DICT_DELETION:
        try:
            names = getattr(module_class, entry)
        except AttributeError:
            continue
        for name, submodule in names.items():
            super(LazyModule, module).__setattr__(name, submodule)


def _reset_lazymodule(module, class_attrs):
    module_class = type(module)
    del module_class.__getattribute__
    del module_class.__setattr__
    try:
        del module_class._LOADING
    except AttributeError:
        pass

    for attr in _CLS_ATTRS:
        try:
            setattr(module_class, attr, class_attrs[attr])
        except KeyError:
            pass

    _reset_lazy_submodule_refs(module)


def _load_module(module):
    module_class = type(module)
    if not issubclass(module_class, LazyModule):
        raise TypeError('Not an instance of LazyModule')
    with _ImportLockContext():
        parent, _, module_name = module.__name__.rpartition('.')
        if not hasattr(module_class, '_lazy_import_error_msgs'):
            return
        module_class._LOADING = True
        try:
            if parent:
                setattr(sys.modules[parent], module_name, module)
            if not hasattr(module_class, '_LOADING'):
                return
            cached_data = _clean_lazymodule(module)
            try:
                reload_module(module)
            except Exception:
                _reset_lazymodule(module, cached_data)
                raise
            else:
                delattr(module_class, '_LOADING')
                _reset_lazy_submodule_refs(module)
        except (AttributeError, ImportError):
            msg = module_class._lazy_import_error_msgs['msg']
            raise_from(ImportError(msg.format(**module_class._lazy_import_error_strings)), None)


def _lazy_module(module_name, error_strings):
    with _ImportLockContext():
        full_module_name = module_name
        full_submodule_name = None
        submodule_name = ''
        while module_name:
            try:
                module = sys.modules[module_name]
                module_name = ''
            except KeyError:
                err_strs = error_strings.copy()
                err_strs.setdefault('module', module_name)

                class _LazyModule(LazyModule):
                    _lazy_import_error_msgs = {'msg': err_strs.pop('msg')}
                    msg_callable = err_strs.pop('msg_callable', None)
                    if msg_callable:
                        _lazy_import_error_msgs['msg_callable'] = msg_callable
                    _lazy_import_error_strings = err_strs
                    _lazy_import_callables = {}
                    _lazy_import_submodules = {}

                    def __repr__(self):
                        return 'Lazily-loaded module %s' % self.__name__

                _LazyModule.__name__ = 'module'
                module = sys.modules[module_name] = _LazyModule(module_name)

            if full_submodule_name:
                submodule = sys.modules[full_submodule_name]
                ModuleType.__setattr__(module, submodule_name, submodule)
                _LazyModule._lazy_import_submodules[submodule_name] = submodule

            full_submodule_name = module_name
            module_name, _, submodule_name = module_name.rpartition('.')

        return sys.modules[full_module_name]


def lazy_module(module_name, level='leaf'):
    module_base_name = module_name.partition('.')[0]
    error_strings = {}
    try:
        caller = sys._getframe(3).f_globals['__name__']
    except AttributeError:
        caller = 'Python'
    error_strings.setdefault('caller', caller)
    error_strings.setdefault('install_name', module_base_name)
    error_strings.setdefault('msg', _MSG)

    module = _lazy_module(module_name, error_strings)

    if level == 'base':
        return sys.modules[module_base_name]
    elif level == 'leaf':
        return module
    else:
        raise ValueError('Must be base or leaf')
