Install the NOMAD client library
================================

We release the NOMAD client library as a Python `distutils <https://docs.python.org/3/library/distutils.html>`_ source distribution.
You can install it the usual way using *pip* (or *conda*).

.. code-block:: sh

    pip install https::/repository.nomad-coe.eu/app/dist/nomad-0.8.0.tar.gz

There are different layers of dependencies that you have to install, in order to use
certain functions of NOMAD. The base install above, will only install the
necessary packages for accessing the NOMAD Archive and use the NOMAD metainfo (see
:ref:`access the archive <access-the-archive-label>`).

Other functions, e.g. using the NOMAD parsers to parse your code output, require
additional dependencies. You can use the ``[extra]`` notation to install these extra
requirements:

.. code-block:: sh

    pip install https::/repository.nomad-coe.eu/app/dist/nomad-0.8.0.tar.gz[parsers]
    pip install https::/repository.nomad-coe.eu/app/dist/nomad-0.8.0.tar.gz[infra]
    pip install https::/repository.nomad-coe.eu/app/dist/nomad-0.8.0.tar.gz[dev]

The various *extras* have the following meaning:

- ``parsing``, everything necessary to run the parsers
- ``infra``, everything to run NOMAD services
- ``dev``, additional tools that are necessary to develop NOMAD
