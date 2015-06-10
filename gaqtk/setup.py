# Copyright (C) 2015 by Per Unneberg
from setuptools import setup, find_packages
import glob
import versioneer

versioneer.VCS = 'git'
versioneer.versionfile_source = 'gaqtk/_version.py'
versioneer.versionfile_build = 'gaqtk/_version.py'
versioneer.tag_prefix = ''  # tags are like 1.2.0
versioneer.parentdir_prefix = 'gaqtk-'  # dirname like 'myproject-1.2.0'

setup(name="gaqtk",
      # version=versioneer.get_version(),
      # cmdclass=versioneer.get_cmdclass(),
      author="Per Unneberg",
      author_email="per.unneberg@scilifelab.se",
      description="Genome assembly quality visualization toolkit",
      license="MIT",
      scripts=glob.glob('scripts/*.py') + glob.glob('bin/*'),
      install_requires=[
          "pyyaml",
          "bokeh",
          "jinja2",
          "nose",
      ],
      test_suite='nose.collector',
      packages=find_packages(exclude=['ez_setup', 'test*']),
      namespace_packages=[
          'gaqtk',
      ],
      package_data={
          'gaqtk': [
              'static/*',
              'templates/*',
          ],
      })
