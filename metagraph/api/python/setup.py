#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

import email.parser
import json
import os
import sys
from setuptools import setup, find_packages

# every file is read beside this script, not from the current directory: the package builds
# from any directory, and from a copy outside the repository (an sdist, a copied folder)
HERE = os.path.dirname(os.path.abspath(__file__))


def _version():
    """The repository's version (package.json at its top) in a checkout; in a copy without
    it, the version its PKG-INFO states (an sdist, or a folder unpacked from one); else
    0+unknown, said on stderr: nothing states a version there."""
    repository = os.path.join(HERE, '..', '..', '..', 'package.json')
    if os.path.exists(repository):
        with open(repository) as f:
            return json.load(f)['version']
    pkg_info = os.path.join(HERE, 'PKG-INFO')
    if os.path.exists(pkg_info):
        with open(pkg_info, encoding='utf-8') as f:
            version = email.parser.Parser().parse(f, headersonly=True).get('Version')
        if version:
            return version
    print('setup.py: neither %s nor %s exists: the version is unknown (0+unknown)'
          % (os.path.normpath(repository), pkg_info), file=sys.stderr)
    return '0+unknown'


with open(os.path.join(HERE, 'README.rst')) as readme_file:
    readme = readme_file.read()

setup_requirements = ['pytest-runner']

with open(os.path.join(HERE, 'requirements.txt')) as f:
    requirements = list(f.readlines())

test_requirements = ['pytest']

version = _version()

setup(
    author="ratschlab",
    author_email='grlab@ratschlab.org',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: 3 :: Only',
    ],
    # metagraph.traverse uses dataclass(slots=True) (3.10)
    python_requires='>=3.10',
    description="Metagraph Toolkit",
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    keywords='MetaGraph',
    name='metagraph-api',
    # metagraph.* too: the traverse subpackage must survive a non-editable install
    packages=find_packages(include=['metagraph', 'metagraph.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/ratschlab/metagraph',
    version=version,
    zip_safe=False,
)
