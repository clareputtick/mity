#!/bin/bash

# update the requirements.txt
pip-compile 
python3 setup.py sdist
python3 setup.py bdist_wheel

version=$(grep version mitylib/_version.py | cut -f2 -d = | sed 's/[\", ]//g')

twine check dist/mitywgs-${version}*

function test {
  twine upload \
    --repository-url https://test.pypi.org/legacy/ \
    --verbose \
    --non-interactive \
    -u drmjc \
    dist/mitywgs-${version}*
}
function public {
  twine upload --verbose --non-interactive -u drmjc dist/mitywgs-${version}*
}
#test
public
