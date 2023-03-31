# Install the Python library

NOMAD provides a Python package called `nomad-lab`.
The package is hosted on [pypi](https://pypi.org/project/nomad-lab/)
and you can install it with *pip* (or conda).

To install the latest stable pypi release, simply use pip:
```sh
pip install nomad-lab
```

!!! attention
    The latest versions might still be considered beta and might not be published to
    pypi. If you require specific new features you might need to install `nomad-lab`
    from our GitLab package registry. To use features of a specific commit or
    branch, consider to [clone and build the project](./develop/setup.md) yourself.


To install the latest release developer releases from our GitLab use:
```sh
pip install nomad-lab --extra-index-url https://gitlab.mpcdf.mpg.de/api/v4/projects/2187/packages/pypi/simple
```

To install an older version of NOMAD (e.g. v0.10.x), you can of use reference
the respective version on pypy:
```sh
pip install nomad-lab==v0.10.12
```

There are different layers of dependencies that you have to install, in order to use certain functions of NOMAD.
The base install above, only installs the necessary packages for
accessing the NOMAD Archive and using the NOMAD metainfo (see access the archive).

Other functions, e.g. using the NOMAD parsers to parse your code output, require additional dependencies.
You can use the [extra] notation to install these extra requirements:

```
pip install nomad-lab[parsing]
pip install nomad-lab[infrastructure]
pip install nomad-lab[dev]
```
The various extras have the following meaning:

- *parsing*, everything necessary to run the parsers

- *infrastructure*, everything to run NOMAD services

- *dev*, additional tools that are necessary to develop NOMAD

