"""
init file for base package.
"""
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="rfmpy",
    version="0.0.1",
    author="Konstantinos Michailos",
    author_email="konstantinos.michailos@gmail.com",
    description="Set of codes for calculating receiver functions and time to depth migration.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    license_files = ('license.txt'),
    packages=['rfmpy.core', 'rfmpy.utils', 'rfmpy.visualisation'],
    zip_safe=False
)
