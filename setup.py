import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="fings",
    version="1.6.3",
    author="Christopher Wardell",
    author_email="github@cpwardell.com",
    description="Filters for Next Generation Sequencing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/cpwardell/FiNGS",
    packages=['fings','fingsexampledata'],
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
