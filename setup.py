from distutils.core import setup

setup(name='cancermuts',
      version='0.1',
      description='Cancer mutation annotation toolkit',
      author='Matteo Tiberti',
      author_email='matteo.tiberti@gmail.com',
      url='https://www.github.com/ELELAB/cancermuts',
      packages=['cancermuts'],
      install_requires=['requests',
                        'bioservices>=1.10.0',
      		        'myvariant',
      		        'pyliftover',
                        'Bio',
                        'bravado',
                        'matplotlib',
                        'pandas',
                        'parse',
                        'urllib3',
                        'future']
     )

