from setuptools import setup, find_packages

setup(
    name="bsv",
    version="0.1",
    packages=find_packages(),
    package_data={'bsv': ['examples/*.py']},
    install_requires=[
        'numpy',
        'requests',
        'matplotlib',  # Add any other dependencies here
    ],
    author="Julie Fabre",
    author_email="julie.mfabre@gmail.com",
    description="A package for brain structure visualization",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Julie-Fabre/brain_street_view",  # if you have a GitHub repository
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)