import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="fings",
    version="1.7.2",
    author="Christopher Wardell",
    author_email="github@cpwardell.com",
    description="Filters for Next Generation Sequencing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/cpwardell/FiNGS",
    packages=['fings'],
    include_package_data=True,
    install_requires=['pyvcf',
                      'pysam',
                      'numpy',
                      'scipy==1.2',
                      'pandas',
                      'joblib',
                      'seaborn',
                      'statsmodels',
                      'editdistance'],
    entry_points={
            "console_scripts": [
            "fings = fings.FiNGS:main"
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='==3.7',
)
