import os
from setuptools import setup

version_py = os.path.join(os.path.dirname(__file__), 'artic', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','').strip()
long_description = """
``artic`` is a pipeline for working with virus sequencing data sequenced with nanopore
"""

HERE = os.path.dirname(__file__)

with open(os.path.join(HERE, "requirements.txt"), "r") as f:
    install_requires = [x.strip() for x in f.readlines()]

setup(
    name="artic",
    version=version,
    install_requires=install_requires,
    requires=['python (>=3.5)'],
    packages=['artic', 'artic.scripts'],
    author="Nick Loman",
    description='A toolset for working with nanopore sequencing data',
    long_description=long_description,
    url="http://poretools.readthedocs.org",
    package_dir={'artic': "artic"},
    package_data={'artic': []},
    zip_safe=False,
    include_package_data=True,
    #scripts=['poretools/scripts/poretools'],
    entry_points={
        'console_scripts': [
            'artic=artic.pipeline:main',
            'align_trim=artic.align_trim:main',
            'margin_cons=artic.margin_cons:main',
            'vcfextract=artic.vcfextract:main',
        ],
    },
    author_email="n.j.loman@bham.ac.uk",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
    )
