#!/usr/bin/env python
"""Setup script to install sospex."""

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

config = {
        'name':'sospex',
        'version':'0.34.beta',
        'description':'SOFIA SPectral EXplorer',
        'long_description':long_description,
        'long_description_type':"text/markdown",
        'author':'Dario Fadda',
        'author_email':'darioflute@gmail.com',
        'url':'https://github.com/darioflute/sospex.git',
        'download_url':'https://github.com/darioflute/sospex',
        'license':'GPLv3+',
        'packages':['sospex'],
        'scripts':['bin/sospex'],
        'install_requires':[
                'python>=3.6',
                'numpy>=1.11',
                'scipy',
                'astropy>=3.0',
                'matplotlib>=3.0.2'
                ],
        'include_package_data':True,
        'package_data':{'sospex':['copyright.txt','icons/*.png','icons/*.gif',
                                  'help/*.html','yellow-stylesheet.css']},
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
