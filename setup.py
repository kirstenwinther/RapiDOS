import setuptools

with open('requirements.txt', 'r') as f:
    requirements = f.readlines()

with open('README.md', 'r') as f:
    readme = f.read()

setuptools.setup(
    name="RapiDOS",
    version="0.1",
    #url="https://github.com/SUNCAT-Center/",
    author="Kirsten Winther",
    author_email="winther@stanford.edu",

    description="Quick (P)DOS plotting and parsing from VASP",
    long_description=readme,
    license='GPL-3.0',

    packages=[
        'rapidos',
    ],
    package_dir={'rapidos': 'rapidos'},
    entry_points={'console_scripts': ['rapidos=rapidos.cli:cli']},
    install_requires=requirements,
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.11',
    ],
)
