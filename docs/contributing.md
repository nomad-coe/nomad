# Contributing

The are some *rules* or better strong *guidelines*

- Use an IDE (e.g. [vscode](https://code.visualstudio.com/)) or otherwise automatically
enforce code [formatting and linting](https://code.visualstudio.com/docs/python/linting).

- There is a style guide to python. Write [pep-8](https://www.python.org/dev/peps/pep-0008/)
compliant python code. An exception is the line cap at 79, which can be broken but keep it 90-ish.

- Test the public API of each submodule (i.e. python file)

- Be [pythonic](https://docs.python-guide.org/writing/style/) and watch
[this](https://www.youtube.com/watch?v=wf-BqAjZb8M).

- Document any *public* API of each submodule (e.g. python file). Public meaning API that
is exposed to other submodules (i.e. other python files).

- Use google [docstrings](http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html).

- Add your docstrings to the sphinx documentation in `docs`. Use .md, follow the example.
Markdown in sphix is supported via [recommonmark](https://recommonmark.readthedocs.io/en/latest/index.html#autostructify)
and [AutoStructify](http://recommonmark.readthedocs.io/en/latest/auto_structify.html)

- The project structure is according to [this](https://docs.python-guide.org/writing/structure/)
guide. Keep it!
