"""
init file for base package.
"""
from setuptools import setup, find_packages

setup(
    name="rfmpy",
    version="0.0.1",
    author="Konstantinos Michailos",
    author_email="konstantinos.michailos@gmail.com",
    description="A small set of codes for calculating receiver functions and time to depth migration.",
    url="",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    license_files = ('license.txt'),
    packages=['rfmpy.core', 'rfmpy.utils', 'rfmpy.visuals'],
    zip_safe=False
)