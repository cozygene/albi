
from distutils.core import setup

setup(
      name='albi',
      version='1.0',
      url = 'https://github.com/cozygene/albi',
      description = '',
      packages=['albi_lib'],
      scripts=['albi_lib/fiesta.py', 'albi_lib/albi.py', 'albi_lib/progress_bar.py'],
    )
