# -*- coding: utf-8 -*-
# Installation script for stl2pov and friends
#
# R.F. Smith <rsmith@xs4all.nl>
# $Date: 2012-06-04 21:17:53 +0200 $

from distutils.core import setup

with open('README.txt') as f:
    ld = f.read()


setup(name='py-stl',
      version='$Revision: 3.1 $'[11:-2],
      license='BSD',
      description='Programs to read and convert STL files.',
      author='Roland Smith', author_email='rsmith@xs4all.nl',
      url='http://rsmith.home.xs4all.nl/software/',
      scripts=['stl2pov.py', 'stl2ps.py', 'stl2pdf.py', 'stlinfo.py'],
      provides='pystl', py_modules=['stl','xform'],
      classifiers=['Development Status :: 5 - Production/Stable',
                   'Environment :: Console',
                   'Intended Audience :: End Users/Desktop',
                   'Intended Audience :: Manufacturing',
                   'License :: OSI Approved :: BSD License',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python :: 2.7',
                   'Topic :: Scientific/Engineering'
                   ],
      long_description = ld
      )
