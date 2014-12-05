#!/usr/bin/env python

from distutils.core import setup

setup(name='waveloc',
	version='0.3.0',
	description='Earthquake location by waveform migration',
	author='Alessia Maggi',
	author_email='alessia.maggi@unistra.fr',
	url='http://github.com/amaggi/waveloc',
	packages=['waveloc', 'waveloc_examples'],
	package_dir = {'waveloc' : 'waveloc',
                   'waveloc_examples' : 'waveloc_examples'},
	scripts=['scripts/grid2hdf5', 'scripts/pyms2sac', 'scripts/make_SDS'],
        requires=[
		'numpy(>=1.6.1)',
		'obspy.core(>=0.7.1)',
		'h5py(>=2.0.0)',
	],
	classifiers=[
		'Environment :: Console',
		'Intended Audience :: Science/Research',
		'License :: OSI Approved :: CeCILL',
		'Programming Language :: Python :: 2.7',
		'Topic :: Scientific/Engineering',
	], 
)
