from setuptools import setup
import os
import re


def read_version():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'scibacr/__init__.py')
    with open(path, 'r') as fh:
        return re.search(r'__version__\s?=\s?[\'"](.+)[\'"]', fh.read()).group(1)


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='scibacr',
      version=read_version(),
      description='',
      long_description=readme(),
      long_description_content_type='text/markdown',
      author='Ariane Mora',
      author_email='ariane.n.mora@gmail.com',
      url='https://github.com/ArianeMora/scibacr',
      license='GPL3',
      project_urls={
          "Bug Tracker": "https://github.com/ArianeMora/scibacr/issues",
          "Documentation": "https://github.com/ArianeMora/scibacr",
          "Source Code": "https://github.com/ArianeMora/scibacr",
      },
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      keywords='util',
      packages=['scibacr'],
      entry_points={
          'console_scripts': [
              'scibacr = scibacr.__main__:main'
          ]
      },
      install_requires=['pandas', 'numpy', 'pysam', 'tqdm', 'h5py', 'sciviso', 'statsmodels',
                        'scipy', 'scivae>=1.0.8', 'wordcloud'],
      python_requires='>=3.8',
      data_files=[("", ["LICENSE"])]
      )