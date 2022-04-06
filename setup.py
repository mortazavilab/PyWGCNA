from setuptools import setup

setup(
    name='PyWGCNA',  # the name of your package
    packages=['PyWGCNA'],  # same as above
    version='v0.5.9',  # version number
    license='MIT',  # license type
    description='PyWGCNA is a Python package designed to do Weighted correlation network analysis (WGCNA)',
    # short description
    author='Narges Rezaie',  # your name
    author_email='nargesrezaie80@gmail.com',  # your email
    url='https://github.com/mortazavilab/PyWGCNA',  # url to your git repo
    download_url='https://github.com/mortazavilab/PyWGCNA/archive/refs/tags/v0.2.4-alpha.tar.gz',  # link to the tar.gz file associated with this release
    keywords=['PyWGCNA', 'WGCNA', 'bulk', 'gene clustering', 'network analysis'],  #
    install_requires=[  # these can also include >, <, == to enforce version compatibility
        'pandas',  # make sure the packages you put here are those NOT included in the
        'numpy',  # base python distro
        'scipy',
        'scikit-learn>=0.24.2',
        'statsmodels',
        'matplotlib',
        'seaborn',
        'biomart',
        'gseapy',
        'pytest',
        'pyvis'
    ],
    classifiers=[  # choose from here: https://pypi.org/classifiers/
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research ',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
    ],
)
