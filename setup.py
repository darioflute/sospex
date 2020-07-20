#!/usr/bin/env python
"""Setup script to install sospex."""

#from distutils.core import setup
import setuptools
from setuptools import setup
import json

with open('sospex/version.json') as fp:
    _info = json.load(fp)

with open("README.md", "r") as fh:
    long_description = fh.read()

config = {
        'name':'sospex',
        'version':_info['version'],
        'description':'SOFIA SPectral EXplorer',
        'long_description':long_description,
        'long_description_type':"text/markdown",
        'author':'Dario Fadda',
        'author_email':'darioflute@gmail.com',
        'url':'https://github.com/darioflute/sospex.git',
        'download_url':'https://github.com/darioflute/sospex',
        'python_requires':'>=3.7',
        'license':'GPLv3+',
        'packages':['sospex'],
        'scripts':['bin/sospex'],
        'include_package_data':True,
        'package_data':{'sospex':['copyright.txt','icons/*.png','icons/*.gif',
                                  'help/*.html','yellow-stylesheet.css',
                                  'data/*.gz','version.json']},
        'classifiers':[
                "Programming Language :: Python :: 3",
                "License :: OSI Approved :: GPLv3+ License",
                "Operating System :: OS Independent",
                "Intended Audience :: Science/Research", 
                "Development Status :: 4 - Beta",
                "Topic :: Scientific/Engineering :: Visualization",
                ],
     }

setup(**config)
