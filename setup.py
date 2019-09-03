import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="idbd_utils",
    version="0.0.1",
    author="IDbyDNA",
    author_email="jmontgomery@idbydna.com",
    description="IDbyDNA utilities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/idbydna/taxonomer",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: None",
        "Operating System :: OS Independent",
    ],
)