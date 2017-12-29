#!/usr/bin/env python

from distutils.core import setup

setup(name='sospex',
      version='0.9.beta',
      description='SOFIA SPectrum EXplorer',
      long_description='The package displays FIFI-LS cubes',
      author='Dario Fadda',
      author_email='darioflute@gmail.com',
      url='https://github.com/darioflute/sospex.git',
      download_url='https://github.com/darioflute/sospex',
      license='GPLv3+',
      packages=['sospex'],
      scripts=['bin/sospex'],
      include_package_data=True,
      package_data={'sospex':['copyright.txt','icons/*.png','icons/*.gif','help/*.html','yellow-stylesheet.css']}
     )
