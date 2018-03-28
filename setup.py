from setuptools import setup, find_packages
import os

setup(
      name="arachis",
      version="0.41.2",
      description="Analyzing ReArrangements with Circular Headed Incomplete Sequences",
      author="Jianjun Jin",
      author_email='jinjianjun@mail.kib.ac.cn',
      url="http://github???",
      license="GNU General Public License, version 3",
      packages=["arachis"],
      scripts=["scripts/run_pypmag.py", "scripts/compare_grimm_files.py", "scripts/map_tree_labels.py"],
      )

os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')