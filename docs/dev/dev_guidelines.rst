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

The CI/CD will run on all refs that do not start with ``dev-``. The CI/CD will
not release or deploy anything automatically, but it can be manually triggered after the
build and test stage completed successfully.


Git/GitLab
^^^^^^^^^^

The ``master`` branch of our repository is *protected*. You must not (even if you have
the rights) commit to it directly. The ``master`` branch references the latest official
release (i.e. what the current NOMAD runs on). The current development is represented by
*version* branches, named ``vx.x.x``. Usually there are two or more of these branched,
representing the development on *minor/bugfix* versions and the next *major* version(s).
Ideally these *version* branches are also not manually push to.

Instead you develop
on *feature* branches. These are branches that are dedicated to implement a single feature.
They are short lived and only exist to implement a single feature.

The lifecycle of a *feature* branch should look like this:

- create the *feature* branch from the last commit on the respective *version* branch that passes CI

- do your work and push until you are satisfied and the CI passes

- create a merge request on GitLab

- discuss the merge request on GitLab

- continue to work (with the open merge request) until all issues from the discussion are resolved

- the maintainer performs the merge and the *feature* branch gets deleted

While working on a feature, there are certain practices that will help us to create
a clean history with coherent commits, where each commit stands on its own.

.. code-block:: sh

  git commit --amend

If you committed something to your own feature branch and then realize by CI that you have
some tiny error in it that you need to fix, try to amend this fix to the last commit.
This will avoid unnecessary tiny commits and foster more coherent single commits. With `amend`
you are basically adding changes to the last commit, i.e. editing the last commit. If
you push, you need to force it ``git push origin feature-branch -f``. So be careful, and
only use this on your own branches.

.. code-block:: sh

  git rebase <version-branch>

Lets assume you work on a bigger feature that takes more time. You might want to merge
the version branch into your feature branch from time to time to get the recent changes.
In these cases, use rebase and not merge. Rebase puts your branch commits in front of the
merged commits instead of creating a new commit with two ancestors. It basically moves the
point where you initially branched away from the version branch to the current position in
the version branch. This will avoid merges, merge commits, and generally leave us with a
more consistent history.  You can also rebase before create a merge request, basically
allowing for no-op merges. Ideally the only real merges that we ever have, are between
version branches.

.. code-block:: sh

  git merge --squash <other-branch>

When you need multiple branches to implement a feature and merge between them, try to
use `squash`. Squashing basically puts all commits of the merged branch into a single commit.
It basically allows you to have many commits and then squash them into one. This is useful
if these commits where just made for synchronization between workstations or due to
unexpected errors in CI/CD, you needed a save point, etc. Again the goal is to have
coherent commits, where each commits makes sense on its own.

Often a feature is also represented by an *issue* on GitLab. Please mention the respective
issues in your commits by adding the issue id at the end of the commit message: `My message. #123`.

We tag releases with ``vX.X.X`` according to the regular semantic versioning practices.
After releasing and tagging the *version* branch is removed. Do not confuse tags with *version* branches.
Remember that tags and branches are both Git references and you can accidentally pull/push/checkout a tag.

The main NOMAD GitLab-project (``nomad-fair``) uses Git-submodules to maintain its
parsers and other dependencies. All these submodules are places in the `/dependencies`
directory. There are helper scripts to install (`dependencies.sh`, see :ref:`setup </setup.html>`) and
commit changes to all submodules (`dependencies-git.sh`). After merging or checking out,
you have to make sure that the modules are updated to not accidentally commit old
submodule commits again. Usually you do the following to check if you really have a
clean working directory.

.. code-block:: sh

  git checkout something-with-changes
  git submodule update
  git status


Terms and Identifiers
---------------------

There are is some terminology consistently used in this documentation and the source
code. Use this terminology for identifiers.

Do not use abbreviations. There are (few) exceptions: ``proc`` (processing); ``exc``, ``e`` (exception);
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


.. _id-reference-label:

Ids
---

Throughout nomad, we use different ids. If something
is called *id*, it is usually a random uuid and has no semantic connection to the entity
it identifies. If something is called a *hash* than it is a hash build based on the
entity it identifies. This means either the whole thing or just some properties of
said entities.

- The most common hashes is the ``calc_hash`` based on mainfile and auxfile contents.
- The ``upload_id`` is a UUID assigned at upload time and never changed afterwards.
- The ``mainfile`` is a path within an upload that points to a main code output file.
  Since, the upload directory structure does not change, this uniquely ids a calc within the upload.
- The ``calc_id`` (internal calculation id) is a hash over the ``mainfile`` and respective
  ``upload_id``. Therefore, each `calc_id` ids a calc on its own.
- We often use pairs of `upload_id/calc_id`, which in many context allow to resolve a calc
  related file on the filesystem without having to ask a database about it.
- The ``pid`` or (``coe_calc_id``) is an sequential interger id.
- Calculation ``handle`` or ``handle_id`` are created based on those ``pid``.
  To create hashes we use :py:func:`nomad.utils.hash`.


NOMAD-coe Dependencies
----------------------

We currently use git submodules to maintain references to NOMAD-coe dependencies.
All dependencies are python packages and installed via pip to your python environement.

This allows us to target (e.g. install) individual commits. More importantly, we can address c
ommit hashes to identify exact parser/normalizer versions. On the downside, common functions
for all dependencies (e.g. the python-common package, or nomad_meta_info) cannot be part
of the nomad-FAIRDI project. In general, it is hard to simultaneously develop nomad-FAIRDI
and NOMAD-coe dependencies.

Another approach is to integrate the NOMAD-coe sources with nomad-FAIRDI. The lacking
availability of individual commit hashes, could be replaces with hashes of source-code
files.

We use the branch ``nomad-fair`` on all dependencies for nomad-FAIRDI specific changes.


Parsers
^^^^^^^

There are several steps to take, to wrap a NOMAD-coe parser into a nomad@FAIRDI parser:

- Implement ``nomadcore.baseclasses.ParserInterface`` or a class with a similar constructutor
  and `parse` method interface.
- Make sure that the meta-info is
  only loaded for each parse instance, not for each parser run.
- Have a root package that bears the parser name, e.g. ``vaspparser``
- The important classes (e.g. the parser interface implementation) in the root module
  (e.g. ``vaspparser/__init__.py``)
- Only use sub-modules were necessary. Try to avoid sub-directories
- Have a test module. Don't go overboard with the test data.
- Make it a pypi-style package, i.e. create ``setup.py`` script.
- The package name should be the parser name, e.g. ``vaspparser``.
- Let the parser logging as it is. We will catch it with a handler installed on the root logger.
  This handler will redirect all legacy log events and put it though the nomad@FAIRDI
  treatment described below.
- Remove all scala code.


Normalizers
^^^^^^^^^^^

We are rewriting all NOMAD-coe normalizers, see :py:mod:`nomad.normalizing`.


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
- We use logging to inform about the state of nomad-FAIRDI, not about user
  behavior, input, data. Do not confuse this when determining the log-level for an event.
  For example, a user providing an invalid upload file, for example, should never be an error.

Please follow the following rules when logging:

- If a logger is not already provided, only use
  :py:func:`nomad.utils.get_logger` to acquire a new logger. Never use the
  build-in logging directly. These logger work like the system loggers, but
  allow you to pass keyword arguments with additional context data. See also
  the `structlog docs <https://structlog.readthedocs.io/en/stable/>`_.
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

Used log keys
^^^^^^^^^^^^^
The following keys are used in the final logs that are piped to Logstash.
Notice that the key name is automatically formed by a separate formatter and
may differ from the one used in the actual log call.

Keys that are autogenerated for all logs:

 - ``@timestamp``: Timestamp for the log
 - ``@version``: Version of the logger
 - ``host``: The host name from which the log originated
 - ``path``: Path of the module from which the log was created
 - ``tags``: Tags for this log
 - ``type``: The `message_type` as set in the LogstashFormatter
 - ``level``: The log level: ``DEBUG``, ``INFO``, ``WARNING``, ``ERROR``
 - ``logger_name``: Name of the logger
 - ``nomad.service``: The service name as configured in ``config.py``
 - ``nomad.release``: The release name as configured in ``config.py``

Keys that are present for events related to processing an entry:

 - ``nomad.upload_id``: The id of the currently processed upload
 - ``nomad.calc_id``: The id of the currently processed entry
 - ``nomad.mainfile``: The mainfile of the currently processed entry

Keys that are present for events related to exceptions:

 - ``exc_info``: Stores the full python exception that was encountered. All
   uncaught exceptions will be stored automatically here.
 - ``digest``: If an exception was raised, the last 256 characters of the message
   are stored automatically into this key. If you wish to search for exceptions
   in Kibana, you will want to use this value as it will be indexed unlike the
   full exception object.
