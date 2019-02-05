Development guidelines
======================

Design principles
-----------------

- simple first, complicated only when necessary
- adopting generic established 3rd party solutions before implementing specific solutions
- only uni directional dependencies between components/modules, no circles
- only one language: Python (except, GUI of course)


Source code & Git repository
----------------------------

Code Rules
^^^^^^^^^^

The are some *rules* or better strong *guidelines* for writing code. The following
applies to all python code (and were applicable, also to JS and other code):

- Use an IDE (e.g. `vscode <https://code.visualstudio.com/>`_ or otherwise automatically
  enforce code (`formatting and linting <https://code.visualstudio.com/docs/python/linting>`_).
  Use ``nomad qa`` before committing. This will run all tests, static type checks, linting, etc.

- There is a style guide to python. Write `pep-8 <https://www.python.org/dev/peps/pep-0008/>`_
  compliant python code. An exception is the line cap at 79, which can be broken but keep it 90-ish.

- Test the public API of each sub-module (i.e. python file)

- Be `pythonic <https://docs.python-guide.org/writing/style/>`_ and watch
  `this <https://www.youtube.com/watch?v=wf-BqAjZb8M>`_.

- Document any *public* API of each sub-module (e.g. python file). Public meaning API that
  is exposed to other sub-modules (i.e. other python files).

- Use google `docstrings <http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`_.

- Add your doc-strings to the sphinx documentation in ``docs``. Use .md, follow the example.
  Markdown in sphix is supported via `recommonmark
  <https://recommonmark.readthedocs.io/en/latest/index.html#autostructify>`_
  and `AutoStructify <http://recommonmark.readthedocs.io/en/latest/auto_structify.html>`_

- The project structure is according to `this guid <https://docs.python-guide.org/writing/structure/>`_.
  Keep it!


CI/CD
^^^^^

These *guidelines* are partially enforced by CI/CD. As part of CI all tests are run on all
branches; further we run a *linter*, *pep8* checker, and *mypy* (static type checker). You can
run ``nomad qa`` to run all these tests and checks before committing.

Only the CI/CD on ``master`` will create new ``*:latest`` images and allow to deploy.


Git/GitLab
^^^^^^^^^^

The ``master`` branch of our repository is *protected*. You must not (even if you have
the rights) commit to it directly. You develop on *feature* branches and commit changes
to the master branch via *merge requests*. After merge feature branches should be removed.

We tag releases with ``vX.X.X`` according to the regular semantic versioning practices.
There might be branches for older versions to facilitate hot-fixes.


Terms and Identifiers
---------------------

There are is some terminology consistently used in this documentation and the source
code. Use this terminology for identifiers.

Do not use abbreviations. There are (few) exceptions: `proc` (processing); `exc`, `e` (exception);
``calc`` (calculation), ``repo`` (repository), ``utils`` (utilities), and ``aux`` (auxiliary).
Other exceptions are ``f`` for file-like streams and ``i`` for index running variables.
Btw., the latter is almost never necessary in python.

Terms:

- upload: A logical unit that comprises one (.zip) file uploaded by a user.
- calculation: A computation in the sense that is was created by an individual run of a CMS code.
- raw file: User uploaded files (e.g. part of the uploaded .zip), usually code input or output.
- upload file/uploaded file: The actual (.zip) file a user uploaded
- mainfile: The mainfile output file of a CMS code run.
- aux file: Additional files the user uploaded within an upload.
- repo entry: Some quantities of a calculation that are used to represent that calculation in the repository.
- archive data: The normalized data of one calculation in nomad's meta-info-based format.

Ids and Hashes
--------------

Throughout nomad, we use different ids and hashes to refer to entities. If something
is called *id*, it is usually a random uuid and has no semantic connection to the entity
it identifies. If something is called a *hash* than it is a hash build based on the
entity it identifies. This means either the whole thing or just some properties of
said entities.

The most common hashes are the *upload_hash* and *calc_hash*. The upload hash is
a hash over an uploaded file, as each upload usually refers to an individual user upload
(usually a .zip file). The calc_hash is a hash over the mainfile path within an upload.
The combination of upload_hash and calc_hash is used to identify calculations. They
allow us to id calculations independently of any random ids that are created during
processing. To create hashes we use :py:func:`nomad.utils.hash`.


NOMAD-coe Dependencies
----------------------

We currently clone and install NOMAD-coe dependencies *"outside"* the nomad-FAIRDI project
(see :py:mod:`nomad.dependencies`). The installed projects become part of the python
environment and all dependencies are used like regular pipy packages and python modules.

This allows us to target (e.g. install) individual commits. In theory, these might
change during runtime, allowing to update parsers or normalizers on a running nomad.
More importantly, we can address commit hashes to identify exact parser/normalizer versions.
On the downside, common functions for all dependencies (e.g. the python-common package,
or nomad_meta_info) cannot be part of the nomad-FAIRDI project. In general, it is hard
to simultaneously develop nomad-FAIRDI and NOMAD-coe dependencies.

Another approach is to integrate the NOMAD-coe sources with nomad-FAIRDI. The lacking
availability of individual commit hashes, could be replaces with hashes of source-code
files.

We use the branch ``nomad-fair`` on all dependencies for nomad-FAIRDI specific changes.


Parsers
^^^^^^^

There are several steps to take, to make a NOMOAD-coe parser fit for nomad-FAIRDI:

- Implement ``nomadcore.baseclasses.ParserInterface``. Make sure that the meta-info is
  only loaded for each parse instance, not for each parser run.
- Have a root package that bears the parser name, e.g. ``vaspparser``
- The important classes (e.g. the parser interface implementation) in the root module
  (e.g. ``vaspparser/__init__.py``)
- Only use sub-modules were necessary. Try to avoid sub-directories
- Have a test module. Don't go overboard with the test data.
- Make it a pypi-style package, i.e. create ``setup.py`` script.
- The package name should be the parser name, e.g. ``vaspparser``.
- The parser should only use the provided logger
  (:py:func:`nomadcore.baseclasses.ParserInterface::setup_logger`).
  This is important for two reasons. First, our logging uses structured logging and
  all entries are tagged data about parser, upload_ids, mainfiles, etc. This is important
  to make errors easily reproduceable. Second, we store all logs on a parser run to
  be available for end users.
- Keep logging sensible (see logging below). Do not log everything. Do not log massive
  amounts of data. Keep in mind what are errors (as in the parser cannot perform its job)
  and what not (the input is faulty).
- Remove all scala code.


Normalizers
^^^^^^^^^^^

There are several steps to take, to make a NOMOAD-coe normalizer fit for nomad-FAIRDI:

- If written in scala, re-write it in python.
- The normalizer should read from the provided backend. In NOMAD-coe normalizers read
  data from provided serialized dictionaries. Don't do that; we do not want to use such
  a normalizer specific interface.
- Do package, module, and logging related changes as you would for a parser.


Logging
-------

There are three important prerequisites to understand about nomad-FAIRDI's logging:

- All log entries are recorded in a central elastic search database. To make this database
  useful, log entries must be sensible in size, frequence, meaning, level, and logger name.
  Therefore, we need to follow some rules when it comes to logging.
- We use an *structured* logging approach. Instead of encoding all kinds of information
  in log messages, we use key-value pairs that provide context to a log *event*. In the
  end all entries are stored as JSON dictionaries with ``@timestamp``, ``level``,
  ``logger_name``, ``event`` plus custom context data. Keep events very short, most
  information goes into the context.
- We use logging to inform us about the state of nomad-FAIRDI, not about user
  behavior, input, data. Do not confuse this when determining the log-level for an event.
  A user providing an invalid upload file, for example, should never be an error.

Please follow the following rules when logging:

- Only use :py:func:`nomad.utils.get_logger` to acquire a logger. Never use the build-in
  logging directly. These logger work like the system loggers, but allow you to
  pass keyword arguments with additional context data. See also the
  `structlog docs <https://structlog.readthedocs.io/en/stable/>`_.
- In many context, a logger is already provided (e.g. api, processing, parser, normalizer).
  This provided logger has already context information bounded. So it is important to
  use those instead of acquiring your own loggers. Have a look for methods called
  ``get_logger`` or attributes called ``logger``.
- Keep events (what usually is called *message*) very short. Examples are: *file uploaded*,
  *extraction failed*, etc.
- Structure the keys for context information. When you analyse logs in ELK, you will
  see that the set of all keys over all log entries can be quit large. Structure your
  keys to make navigation easier. Use keys like ``nomad.proc.parser_version`` instead of
  ``parser_version``. Use module names as prefixes.
- Don't log everything. Try to anticipate, how you would use the logs in case of bugs,
  error scenarios, etc.
- Don't log sensitive data.
- Think before logging data (especially dicts, list, numpy arrays, etc.).
- Logs should not be abused as a *printf*-style debugging tool.
