from setuptools import setup

setup(name="dynamic_analysis",
      packages=['core'],
      package_data={'core': ['_core.so']})
