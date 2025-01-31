from setuptools import setup, find_packages
import re


def get_version(verbose=1):
    """ Extract version information from source code """

    try:
        with open('mobilityapp/version.py', 'r') as f:
            ln = f.readline()
            # print(ln)
            m = re.search('.* ''(.*)''', ln)
            version = (m.group(1)).strip('\'')
    except Exception as E:
        print(E)
        version = 'none'
    if verbose:
        print('get_version: %s' % version)
    return version


def readme():
    with open('README.md') as f:
        return f.read()

setup(name='mobilityapp',
    version=get_version(),

    maintainer='Damon Carrad',
    maintainer_email='damonc@dtu.dk',
    description='Python-based tool to analyse conductance vs gate voltage data '
              'and extract mobility and density.',
    long_description=readme(),
    url='https://github.com/djcarrad/MobilityApp',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering'
    ],
    license='MIT',

    packages=find_packages(),

    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'PyQt5',
        'ipykernel',
        'lmfit'
    ],
)

version_template = '''
*****
***** package {0} must be at least version {1}.
***** Please upgrade it (pip install -U {0} or conda install {0})
***** in order to use {2}
*****
'''

missing_template = '''
*****
***** package {0} not found
***** Please install it (pip install {0} or conda install {0})
***** in order to use {1}
*****
'''

valueerror_template = '''
*****
***** package {0} version not understood
***** Please make sure the installed version ({1})
***** is compatible with the minimum required version ({2})
***** in order to use {3}
*****
'''

othererror_template = '''
*****
***** could not import package {0}. Please try importing it from 
***** the commandline to diagnose the issue.
*****
'''