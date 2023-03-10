import os
import re
from setuptools import setup, find_packages


def rel_file(*args):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), *args)


def read_from(filename):
    with open(filename) as fp:
        return fp.read()


def get_version():
    data = read_from(rel_file('thesis', '__init__.py'))
    return re.search(r"__version__ = '([^']+)'", data).group(1)


def get_requirements():
    data = read_from(rel_file('requirements.txt'))
    lines = map(lambda s: s.strip(), data.splitlines())
    return list(filter(None, lines))


setup(
    name='thesis',
    version=get_version(),
    description="This library implements all the useful code for Pierre Le Barbenchon's thesis",
    long_description=read_from('README.md'),
    author='Pierre Le Barbenchon',
    packages=find_packages(),
    install_requires=get_requirements(),
)
