#!/bin/bash

# update the requirements.txt
pip-compile 
python3 setup.py sdist
python3 setup.py bdist_wheel

version=$(grep version mitylib/_version.py | cut -f2 -d = | sed 's/[\", ]//g')
function test {
  twine upload --verbose --repository-url https://test.pypi.org/legacy/ dist/mitywgs-${version}* -u drmjc
}
function public {
  twine upload --verbose dist/mitywgs-${version}* -u drmjc
}
test
