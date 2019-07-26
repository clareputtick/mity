#!/bin/bash

# update the requirements.txt
pip-compile 
python setup.py sdist
python setup.py bdist_wheel

version=$(grep version mitylib/_version.py | cut -f2 -d = | sed 's/[\", ]//g')
twine upload --verbose --repository-url https://test.pypi.org/legacy/ dist/mity-${version}* -u drmjc
