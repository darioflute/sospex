#!/usr/bin/env python
"""Setup script to install sospex."""

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

config = {
      'name':'sospex',
      'version':'0.33.beta',
      'description':'SOFIA SPectral EXplorer',
      'long_description':'The package displays FIFI-LS and GREAT cubes',
      'author':'Dario Fadda',
      'author_email':'darioflute@gmail.com',
      'url':'https://github.com/darioflute/sospex.git',
      'download_url':'https://github.com/darioflute/sospex',
      'license':'GPLv3+',
      'packages':['sospex'],
      'scripts':['bin/sospex'],
      'include_package_data':True,
      'package_data':{'sospex':['copyright.txt','icons/*.png','icons/*.gif',
                                'help/*.html','yellow-stylesheet.css']}
     }

setup(**config)
