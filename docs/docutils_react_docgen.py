# =====================
# The MIT License (MIT)
# =====================

# Copyright (c) 2015 Paul Wexler

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

"""
docutils_react_docgen
=====================

docutils extension for documenting React modules.
Requires react-docgen
"""

from docutils import statemachine
from docutils.parsers import rst
import json
import os
import re
import subprocess

REACT_DOCGEN = 'react-docgen' # **Deprecated** Use SETTINGS['react_docgen']

MODULE_DESCRIPTION_MISSING = 'Module doc string is missing!'

MODULE_PROP_DESCRIPTION_MISSING = 'Property doc string is missing!'

MODULE_UNDERLINE_CHARACTER = '-'

TAB_SIZE = 4

DEFAULT_OPTIONS = {
    'exclude':                          '',
    'include':                          '',
    'module_description_missing':       MODULE_DESCRIPTION_MISSING,
    'module_prop_description_missing':  MODULE_PROP_DESCRIPTION_MISSING,
    'module_underline_character':       MODULE_UNDERLINE_CHARACTER,
    'path_index':                       0,
    'show_prop_type':                   False,
    'src':                              '',
    'tab_size':                         TAB_SIZE,
    'use_commonjs_module_name':         True,
}

SETTINGS = {
        'project_base': None,           # absolute path to project base.
        'react_docgen': 'react-docgen', # react-docgen command to run.
        'rst_output': None,             # output filename or no rst output.
        }


def find_package(dirname):
    """find commonjs package for given directory.
    Starts from `dirname` and recurses up the directory tree
    looking for bower.json or package.json.

    Returns a tuple (dirname, package)

    dirname
        The directory the .json file was found in.

    package
        A dict loaded from the .json file.
        Its keys are the module filenames.
    """
    if dirname:
        bower_json = os.path.join(dirname, 'bower.json')
        if os.path.exists(bower_json):
            with open(bower_json, 'r') as f:
                return dirname, json.load(f)
        package_json = os.path.join(dirname, 'package.json')
        if os.path.exists(package_json):
            with open(package_json, 'r') as f:
                return dirname, json.load(f)
        next_dirname = os.path.dirname(dirname)
        if next_dirname != dirname:
            return find_package(next_dirname)
    return None, None


def get_dirname(doc_dict, options):
    return (os.path.dirname(list(doc_dict.keys())[0])
            if options['use_commonjs_module_name'] and doc_dict
            else '')

def react_docgen(args, react_docgen=REACT_DOCGEN):
    """Execute `react-docgen` with the given arguments.

    `args` is a string which may contain spaces.

    `react_docgen` is also a string which may contain spaces.

    Returns the output of `react-docgen` as a dict
    whose keys are module filenames (strings),
    and whose values are module metadata (dicts).

    WARNING
        The default for react_docgen always evaluates to its initial value.
    """
    cmd = react_docgen.split() + args.split()
    return json.loads(subprocess.check_output(cmd, stderr=subprocess.PIPE))

def react_doc_to_rst(doc_dict, options, formatter_class, args=''):
    """ Convert `doc_dict`, the react-docgen output dict
    to a string of ReStructuredText,
    according to the `options` and using the `formatter_class`

    `args` is the string of arguments passed to the directive.
    The default is ''.  Only required when using absolute addressing
    so the Formatter can recover the path_argument.
    """
    dirname = get_dirname(doc_dict, options)
    formatter = formatter_class(options, dirname, args=args)
    return formatter.run(doc_dict)


def run_react_docgen(args, options=DEFAULT_OPTIONS):
    """ Execute `SETTINGS['react_docgen']` with the given args.

    `args` is a string which may contain spaces.

    `SETTINGS['react_docgen']` is also a string which may contain spaces.

    `options` is a dict of directive options.

    The command output is expected to be a JSON blob representing
    a dict whose keys are the module filenames (strings),
    and whose values are the module metadata (dicts).

    However, the blob is simply converted into a python object and returned.

    Implements the `project_base` setting and the `path_index` option.
    """
    arg_list = args.split()
    project_base = SETTINGS['project_base']
    if project_base != None:
        path_index = options['path_index']
        path_argument = arg_list[path_index]
        if not path_argument.startswith(os.path.sep):
            arg_list[path_index] = os.path.abspath(os.path.join(
                    project_base,
                    path_argument))
    cmd = SETTINGS['react_docgen'].split() + arg_list
    return json.loads(subprocess.check_output(cmd, stderr=subprocess.PIPE))

class Formatter(object):
    """ Formatter(options, dirname).run(doc_dict) returns a string.

    options
        a dict of options.

    dirname
        the directory to search for the CommonJS package
        if the use_commonjs_module_name option is True

    doc_dict
        a dict of react-docgen module metadata

    args
        Default is ''.  The string of arguments passed to the directive.
        Required when using absolute addressing.
    """

    def __init__(self, options, dirname, args=''):
        self.options = options
        self.tab = ' ' * self.options['tab_size']
        package_dirname, package = find_package(dirname)
        if package_dirname:
            self.package_dirname_len = len(package_dirname)
            self.package_name = package['name']
        else:
            self.package_dirname_len = 0
        self.args = args
        self._compile_filters()

    def _compile_filters(self):
        include = self.options['include']
        self.include = re.compile(include) if include else None
        exclude = self.options['exclude']
        self.exclude = re.compile(exclude) if exclude else None

    def _filter(self, filename, module_blob):
        """returns True/False to include/exclude the given module
        from the output.
        """
        description = module_blob.get('description', '')
        return ((not self.include or self.include.search(description))
                and
                (not self.exclude or not self.exclude.search(description)))

    def _get_module_name(self, filename):
        if self.package_dirname_len:
            module_name = '%s%s' % (
                    self.package_name,
                    filename[self.package_dirname_len:])
            if module_name.endswith('.js'):
                module_name = module_name[:-3]
        else:
            module_name = filename
        return module_name

    def _get_object_name(self, obj):
        if 'value' in obj:
            value = obj['value']
            if isinstance(value, str):
                return value
            elif isinstance(value, list):
                return '[%s]' % ', '.join(
                        self._get_object_name(item) for item in value)
            else:
                return str(value)
        elif 'name' in obj:
            return obj['name']
        else:
            # if this happens show the obj instead of raising an error.
            return str(obj)

    def _make_definition(self, term, term_definition):
        definition = '\n'.join((self.tab + line
                for line in term_definition.split('\n'))
                if term_definition else [self.tab])
        return term + '\n' + definition

    def _make_emphasis(self, text, style):
        s = ''
        s += style + text + style
        return s

    def _make_heading(self, text, underline_char):
        s = ''
        s += text + '\n'
        s += underline_char * len(text) + '\n\n'
        return s

    def _make_module(self, filename, module_blob):
        s = ''
        s += self._make_module_header(filename)
        s += self._make_module_description(module_blob)
        s += self._make_module_props(module_blob)
        return s

    def _make_module_description(self, module_blob):
        s = ''
        description = module_blob.get('description', '')
        s += description if description else self.options[
                'module_description_missing']
        s += '\n\n'
        return s

    def _make_module_header(self, filename):
        module_name = self._get_module_name(filename)
        s = ''
        s += self._make_heading(
                module_name,
                self.options['module_underline_character'])
        s += self._make_src_link(filename)
        return s

    def _make_module_prop_name(self, name, prop):
        args = []
        if self.options['show_prop_type'] and prop.get('type'):
            args.append(self._get_object_name(prop['type']))
        if prop.get('required'):
            args.append('required')
        if prop.get('defaultValue'):
            args.append('default = ``%s``' % prop['defaultValue']['value'])
        if args:
            return '%s (%s)' % (name, ', '.join(args))
        else:
            return name

    def _make_module_prop(self, name, prop):
        return self._make_definition(
                self._make_module_prop_name(name, prop),
                self._make_module_prop_description(prop)) + '\n\n'

    def _make_module_prop_description(self, prop):
        return prop.get(
                'description',
                self.options['module_prop_description_missing'])

    def _make_module_props(self, module_blob):
        s = ''
        props = module_blob.get('props', {})
        for key in sorted(props.keys()):
            s += self._make_module_prop(key, props[key])
        return s

    def _make_src_link(self, filename):
        s = ''
        if self.options['src']:
            if self.args:
                path = self.args.split()[self.options['path_index']]
                module_name = filename[filename.index(path):]
            else:
                module_name = filename
            link = '%s/%s' % (
                    self.options['src'],
                    module_name)
            s += '`%s`_' % module_name
            s += '\n\n'
            s += '.. _`%s`: %s' % (
                    module_name,
                    link)
            s += '\n\n'
        return s

    def run(self, doc_dict):
        return ''.join(self._make_module(k, d)
                for k, d in sorted(doc_dict.items())
                if self._filter(k, d))


class ReactDocgen(rst.Directive):
    """ Docutils Directive which calls the react-docgen executable.
    """

    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True
    option_spec = {'src': rst.directives.unchanged}
    option_spec.update(
            {k.lower(): rst.directives.unchanged
                    for k in DEFAULT_OPTIONS.keys()})
    has_content = False
    formatter_class = Formatter

    def run(self):
        args = self.arguments[0]
        options = {}
        options.update(DEFAULT_OPTIONS)
        options.update(self.options)
        doc_dict = run_react_docgen(args, options=options)
        rst = react_doc_to_rst(
                doc_dict,
                options,
                self.formatter_class,
                args=args)
        if SETTINGS['rst_output']:
            with open(SETTINGS['rst_output'], 'w') as fo:
                fo.write(rst)
        tab_size = options['tab_size']
        include_lines = statemachine.string2lines(
                rst,
                tab_size,
                convert_whitespace=True)
        self.state_machine.insert_input(include_lines, '')
        return []

rst.directives.register_directive('reactdocgen', ReactDocgen)

