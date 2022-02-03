#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

import json
import os
from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

setup_requirements = ['pytest-runner']

with open('requirements.txt') as f:
    requirements = list(f.readlines())

test_requirements = ['pytest']

with open(os.path.abspath(os.path.dirname(__file__)) + '/../../../package.json') as f:
    version = json.load(f)['version']

setup(
    author="ratschlab",
    author_email='grlab@ratschlab.org',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3 :: Only',
    ],
    description="Metagraph Toolkit",
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    keywords='MetaGraph',
    name='metagraph-api',
    packages=find_packages(include=['metagraph']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/ratschlab/metagraph',
    version=version,
    zip_safe=False,
)
