import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mity",
    version="0.0.1",
    author="Clare Puttick",
    author_email="clare.puttick@gmail.com",
    description="A sensitive Mitochondrial variant detection pipeline from WGS data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KCCG/mity",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)