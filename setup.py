from setuptools import setup

setup(
    name='PyWGCNA',  # the name of your package
    packages=['PyWGCNA'],  # same as above
    version='1.15.4',  # version number
    license='MIT',  # license type
    description='PyWGCNA is a Python package designed to do Weighted correlation network analysis (WGCNA)',
    # short description
    author='Narges Rezaie',  # your name
    author_email='nargesrezaie80@gmail.com',  # your email
    url='https://github.com/mortazavilab/PyWGCNA',  # url to your git repo
    download_url='https://github.com/mortazavilab/PyWGCNA/archive/refs/tags/V1.10.tar.gz',  # link to the tar.gz file associated with this release
    keywords=['PyWGCNA', 'WGCNA', 'bulk', 'gene clustering', 'network analysis'],  #
    install_requires=[  # these can also include >, <, == to enforce version compatibility
        'pandas>=1.3.5',  # make sure the packages you put here are those NOT included in the
        'numpy>=1.21.2',  # base python distribution
        'scipy>=1.9.3',
        'scikit-learn>=1.1.3',
        'statsmodels>=0.12.2',
        'matplotlib>=3.4.2',
        'seaborn>=0.11.1',
        'biomart>=0.9.2',
        'gseapy>=0.14.5',
        'pytest',
        'pyvis>=0.1.9',
        'setuptools>=59.8.0',
        'biomart>=0.9.2',
        'reactome2py>=3.0.0',
        'anndata>=0.8.0',
        'requests>=2.28.1',
    ],
    classifiers=[  # choose from here: https://pypi.org/classifiers/
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research ',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
    ],
)
