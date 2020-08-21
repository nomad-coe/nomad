Using the NOMAD parsers
***********************

To use the NOMAD parsers from the command line, you can use the ``parse`` command. The
parse command will automatically *match* the right parser to your code output file and
run the parser. There are two output formats, ``--show-metadata`` (a JSON representation
of the repository metadata), ``--show-archive`` (a JSON representation of the archive data).

.. code-block:: sh

    nomad parser --show-archive <path-to-your-mainfile-code-output-file>

You can also use the NOMAD parsers from within Python. This will give you the parse
results as metainfo objects to conveniently analyse the results in Python. See :ref:`metainfo <metainfo-label>`
for more details on how to use the metainfo in Python.

.. literalinclude:: ../../examples/parse.py
    :language: python
