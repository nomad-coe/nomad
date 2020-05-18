Install the NOMAD client library
================================

We release the NOMAD client library as a Python `distutils <https://docs.python.org/3/library/distutils.html>`_ source distribution.
You can install it the usual way using *pip* (or *conda*).

.. parsed-literal::

    pip install nomad --extra-index-url |pypi_url|

There are different layers of dependencies that you have to install, in order to use
certain functions of NOMAD. The base install above, will only install the
necessary packages for accessing the NOMAD Archive and use the NOMAD metainfo (see
:ref:`access the archive <access-the-archive-label>`).

Other functions, e.g. using the NOMAD parsers to parse your code output, require
additional dependencies. You can use the ``[extra]`` notation to install these extra
requirements:

.. parsed-literal::

    pip install nomad[parsing] --extra-index-url |pypi_url|
    pip install nomad[infrastructure] --extra-index-url |pypi_url|
    pip install nomad[dev] --extra-index-url |pypi_url|

The various *extras* have the following meaning:

- ``parsing``, everything necessary to run the parsers
- ``infrastructure``, everything to run NOMAD services
- ``dev``, additional tools that are necessary to develop NOMAD
