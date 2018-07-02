from distutils.core import setup

setup(name='cancermuts',
      version='0.1',
      description='Cancer mutation annotation toolkit',
      author='Matteo Tiberti',
      author_email='matteo.tiberti@gmail.com',
      url='',
      packages=['cancermuts'],
      requires=['csv',
      			'requests',
      			'bioservices',
      			'myvariant',
      			'pyliftover',]
     )

