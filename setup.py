from setuptools import setup
try:  # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError:  # for pip <= 9.0.3
    from pip.req import parse_requirements

# parse_requirements() returns generator of pip.req.InstallRequirement objects
install_reqs = parse_requirements('requirements.txt', session='hack')
# reqs is a list of requirement
# e.g. ['django==1.5.1', 'mezzanine==1.4.6']
reqs = [str(ir.req) for ir in install_reqs if 'sphinxcontrib.httpdomain' not in str(ir.req)]

setup(
    name='nomad',
    version='0.5.0',
    description='The nomad@FAIRDI infrastructure python package',
    py_modules=['nomad'],
    install_requires=reqs,
    entry_points='''
        [console_scripts]
        nomad=nomad.client:cli
        admin=nomad.admin:cli
    ''')
