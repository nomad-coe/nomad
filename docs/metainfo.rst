.. _metainfo-label:

Metainfo
========

.. automodule:: nomad.metainfo


Accessing the Metainfo
----------------------

Above you learned what the metainfo is and how to create metainfo definitions and work
with metainfo data in Python. But how do you get access to the existing metainfo definitions
within NOMAD? We call the complete set of all metainfo definitions the *NOMAD Metainfo*.

This *NOMAD Metainfo* comprises definitions from various packages defined by all the
parsers and converters (and respective code outputs and formats) that NOMAD supports. In
addition there are *common* packages that contain definitions that might be relevant to
different kinds of archive data.

Python
______

In the NOMAD source-code all metainfo definitions are materialized as Python source files
that contain the definitions in the format described above. If you have installed the
NOMAD Python package (see :ref:`install-client`), you can simply import the respective
Python modules:

.. code-block:: python

    from nomad.datamodel.metainfo.public import m_package
    print(m_package.m_to_json(indent=2))

    from nomad.datamodel.metainfo.public import section_run
    my_run = section_run()

API
___

In addition, a JSON version of the NOMAD Metainfo is available through our API via the
``metainfo`` endpoint.
You can get :api:`one giant JSON with all definitions <metainfo/>`, or you
can access the metainfo for specific packages, e.g. the :api:`VASP metainfo <metainfo/vasp.json>`. The
returned JSON will also contain all packages that the requested package depends on.

Legacy metainfo version
_______________________

There are no metainfo files anymore. The old ``*.nomadmetainfo.json`` files are no
longer maintained, as the Python definitions in each parser/converter implementation are
now the normative artifact for the NOMAD Metainfo.

To get the NOMAD Metainfo in the format of the old NOMAD CoE project, you can use
the ``metainfo/legacy`` endpoint; e.g. the :api:`VASP legacy metainfo <metainfo/legacy/vasp.nomadmetainfo.json>`.
