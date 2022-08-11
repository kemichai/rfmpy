"""
init file for base package.
"""
from setuptools import setup, find_packages
import os


# Check if we're on the RTD, don't require dependencies if so because we aren't
# actually running the code
READ_THE_DOCS = os.environ.get("READTHEDOCS", None) == "True"

# RTD doesn't require dependencies
if READ_THE_DOCS:
    install_requires = []
else:
    install_requires=[
            # Install requires has been shifted to 'requirements.txt' to avoid
            # Using Pip for package installation. Conda install preferred
            "obspy>=1.2",
            "pandas==1.1.0",
            "matplotlib<3.3",
            "numpy"
            "pyproj>=2.5",
            "scipy",
            ]

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="rfmpy",
    version="0.1.0",
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
