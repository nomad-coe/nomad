## Installing

```
pip install -r requirements.txt
pip install -e .
```

## Contributing

### Code quality

- There is a style guide to python. Write [pep-8](https://www.python.org/dev/peps/pep-0008/) compliant python code. An exception is the line cap at 79, which can be broken but keep it 90-ish.
- Use docstrings and document any *public* API of each submodule (e.g. python file). Public meaning API that is exposed to other submodules (i.e. other python files). Keep docstrings compliant to [pep-257](https://www.python.org/dev/peps/pep-0257/).

### Tests

To run the test execute this from the project root:
```

```

### Generate the documentation

Install [pdoc](https://pypi.org/project/pdoc/) if necessary, and run this from the project root:
```
pdoc --html --html-dir=./build/docs --overwrite nomad
```

Do serve the current docs via HTTP during development, run
```
pdoc --http --http-port 8080
```
and open [http://localhost:8080/nomad](http://localhost:8080/nomad).

