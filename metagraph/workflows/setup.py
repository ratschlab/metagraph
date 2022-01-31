#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

import json
from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

setup_requirements = ['pytest-runner']

with open('requirements.txt') as f:
    requirements = list(f.readlines())

test_requirements = ['pytest']

with open('../../package.json') as f:
    version = json.load(f)['version']

setup(
    author="Marc Zimmermann",
    author_email='marc.zimmermann@inf.ethz.ch',
    maintainer="Mikhail Karasikov",
    maintainer_email='mikhaika@inf.ethz.ch',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Metagraph workflows",
    entry_points={
        'console_scripts': [
            'metagraph-workflows=metagraph_workflows.cli:main'
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    keywords='metagraph_workflows',
    name='metagraph_workflows',
    packages=find_packages(include=['metagraph_workflows']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/ratschlab/metagraph',
    version=version,
    zip_safe=False,
)
