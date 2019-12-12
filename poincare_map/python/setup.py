from setuptools import setup

setup(name="dynamic_analysis",
      packages=['dynamic_analysis', 'dynamic_analysis.core'],
      package_data={'dynamic_analysis.core': ['_core.so']})
