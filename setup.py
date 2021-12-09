from setuptools import setup, find_packages
from os import path


this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), 'r', encoding='utf-8') as f:
    long_description = f.read()

DESCRIPTION = 'Tools to download, parse and write FASTA out of GenBank files',
LONG_DESCRIPTION = long_description,
LONG_DESCRIPTION_CONTENT_TYPE = 'text/markdown'
NAME = 'genbankpy'
AUTHOR = "Semidán Robaina Estévez"
AUTHOR_EMAIL = "srobaina@ull.edu.es"
MAINTAINER = "Semidán Robaina Estévez"
MAINTAINER_EMAIL = "srobaina@gmail.com"
DOWNLOAD_URL = 'http://github.com/robaina/GenBankpy'
LICENSE = 'Creative Commons Attribution 4.0 International'
VERSION = '0.0.1'

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESCRIPTION_CONTENT_TYPE,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=DOWNLOAD_URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      packages=find_packages(),
      install_requires=['ncbi-acc-download', 'biopython'],
    #   entry_points ={
    #         'console_scripts': [
    #             'genbankpy = genbankpy.cli:main'
    #         ]
    #     }
      )