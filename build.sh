#!/bin/bash
python setup.py sdist
python setup.py bdist_wheel
twine upload --verbose --repository-url https://test.pypi.org/legacy/ dist/*
