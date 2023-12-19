# How to install nomad-lab

We provide a Python package called `nomad-lab`. The package can be used to run
certain NOMAD features within local Python programming environments. It includes
the NOMAD parsers and normalizers, or convenience functions to query the processed data on NOMAD.

Released version of the package are hosted on [pypi](https://pypi.org/project/nomad-lab/){:target="_blank"}
and you can install it with *pip* (or conda).

To install the newest pypi release, simply use pip:
```
pip install nomad-lab
```

!!! warning "Attention"
    The latest develop versions might still be considered beta and might not be published to
    pypi. If you require specific new features you might need to install `nomad-lab`
    from our GitLab package registry. To use features of a specific commit or
    branch, consider to [clone and build the project](../develop/setup.md) yourself.


To install the latest release developer releases from our GitLab use:
```
pip install nomad-lab --extra-index-url https://gitlab.mpcdf.mpg.de/api/v4/projects/2187/packages/pypi/simple
```

To install an older version of NOMAD (e.g. v0.10.x), you can of use reference
the respective version on pypy:
```
pip install nomad-lab==1.0.10
```

Certain functionality might require more dependencies. The basic install above,
installs the dependencies for accessing the NOMAD Archive or running most of the NOMAD
parsers.

Other functions, e.g. running the NOMAD infrastructure, require additional dependencies.
You can use the `[extra]` notation to install these extra requirements:

```
pip install nomad-lab[parsing]
pip install nomad-lab[infrastructure]
pip install nomad-lab[dev]
```
The various extras have the following meaning:

- *parsing*, run all parsers, incl. parsers based on HDF5, netCDF, or asr
- *infrastructure*, everything to run NOMAD services, see also
[Oasis documentation](../oasis/install.md#base-linux-without-docker)
- *dev*, necessary to run development and build tools, e.g. pytest, pylint, mypy

