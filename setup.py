from distutils.core import setup

setup(name='alertmap',
      version='0.1dev',
      description='AlertMap Earthquake Early Warning Simulator',
      author='Mike Hearne',
      author_email='mhearne@usgs.gov',
      url='',
      scripts = ['alertmap.py'],
      package_data = {'ttimes.csv']},
)
