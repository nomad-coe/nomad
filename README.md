This project tries and test approaches that might lead to an improved architecture for NOMAD (XT).

## Generate the docs (and continue there)

First, clone this repo:
```
git clone git@gitlab.mpcdf.mpg.de:mscheidg/nomad-xt.git
cd nomad-xt
```

Second, create and source your own virtual python environment:
```
pip install virtualenv
virtualenv -p `which phyton3` .pyenv
source .pyenv/bin/activate
```

Third, install the documentation system [sphinx](http://www.sphinx-doc.org/en/master/index.html):
```
pip install sphinx
pip install recommonmark
```

Forth, generate the documentation:
```
cd docs
make html
```

Conintue with reading the documentation for further setup and contribution guidelines:
```
cd .build/html
python -m http.server 8888
```
Open [http://localhost:8888/html/setup.html](http://localhost:8888/html/setup.html) in
your browser.
