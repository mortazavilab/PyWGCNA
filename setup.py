from setuptools import setup

setup(
    name='PyWGCNA',  # the name of your package
    packages=['PyWGCNA'],  # same as above
    version='2.0.0',  # version number
    license='MIT',  # license type
    description='PyWGCNA is a Python package designed to do Weighted correlation network analysis (WGCNA)',
    # short description
    author='Narges Rezaie',  # your name
    author_email='nargesrezaie80@gmail.com',  # your email
    url='https://github.com/mortazavilab/PyWGCNA',  # url to your git repo
    download_url='https://github.com/mortazavilab/PyWGCNA/archive/refs/tags/v2.0.0.zip',  # link to the tar.gz file associated with this release
    keywords=['PyWGCNA', 'WGCNA', 'bulk', 'gene clustering', 'network analysis'],  #
    install_requires=[  # these can also include >, <, == to enforce version compatibility
        'pandas>=2.1.0',  # make sure the packages you put here are those NOT included in the
        'numpy>=1.24.0',  # base python distribution
        'scipy>=1.9.1',
        'scikit-learn>=1.2.2',
        'statsmodels>=0.14.0',
        'matplotlib>=3.5.2',
        'seaborn>=0.11.2',
        'biomart>=0.9.2',
        'gseapy>=1.0.1',
        'pyvis==0.3.1',
        'setuptools>=67.4.0',
        'biomart>=0.9.2',
        'reactome2py>=3.0.0',
        'anndata>=0.8.0',
        'requests>=2.28.1',
        'networkx>=2.8.4',
        'rsrc>=0.1.3',
        'psutil>=5.9.0',
        'requests>=2.28.1',
    ],
    classifiers=[  # choose from here: https://pypi.org/classifiers/
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research ',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.10',
    ],
)