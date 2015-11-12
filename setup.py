
from distutils.core import setup

setup(
      name='albi',
      version='1.0',
      url = 'https://github.com/cozygene/albi',
      description = '',
      packages=['albi_lib'],
      scripts=['albi.py'],
      install_requires = ['numpy'],
      dependency_links = ['https://pypi.python.org/pypi/numpy#downloads']#TODO Add
    )
