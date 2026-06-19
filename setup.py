import re
from pathlib import Path

from setuptools import setup


def read_version():
    init_file = Path(__file__).parent / 'cancermuts' / '__init__.py'
    init_text = init_file.read_text()
    version_match = re.search(r'^__version__ = [\'"]([^\'"]+)[\'"]', init_text, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError('Unable to find package version.')

setup(name='cancermuts',
      version=read_version(),
      description='Cancer mutation annotation toolkit',
      author='Matteo Tiberti',
      author_email='matteo.tiberti@gmail.com',
      url='https://www.github.com/ELELAB/cancermuts',
      packages=['cancermuts'],
      install_requires=['requests',
                        'bioservices>=1.10.0',
                        'gget',
                        'biothings_client',
                        'pyliftover',
                        'Bio',
                        'bravado',
                        'matplotlib',
                        'pandas',
                        'parse',
                        'urllib3',
                        'future']
     )
