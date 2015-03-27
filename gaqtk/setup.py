# Copyright (C) 2015 by Per Unneberg
from setuptools import setup, find_packages
import os
import glob

setup(name = "gaqtk",
      version = "0.1.0",
      author = "Per Unneberg",
      author_email = "per.unneberg@scilifelab.se",
      description = "Genome assembly quality visualization toolkit",
      license = "MIT",
      scripts = glob.glob('scripts/*.py'),
      install_requires = [
          "pyyaml",
          "mpld3",
          "jinja2",
          "matplotlib",
          "nose",
      ],
      test_suite = 'nose.collector',
      packages=find_packages(exclude=['ez_setup', 'test*']),
      namespace_packages = [
          'gaqtk',
      ],
      package_data = {
        'gaqtk' : [
            'static/*',
            ],
        },
  )

os.system("git rev-parse --short --verify HEAD > ~/.viz_version")