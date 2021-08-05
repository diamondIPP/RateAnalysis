import setuptools

with open('README.md') as f:
    long_description = ''.join(f.readlines())


setuptools.setup(
    name='RateAnalysis',
    version='3.0',
    packages=setuptools.find_packages(),
    include_package_data=True,
    description='Rate Analysis',
    long_description=long_description,
    author='Michael Reichmann',
    author_email='micha.reichmann@gmail.com',
    url='https://github.com/diamondIPP/RateAnalysis',

    install_requires=[
        'ipython',
        'scipy',
        'uncertainties',
        'gtts',
        'pytz',
        'termcolor',
        'progressbar',
        'h5py',
        'screeninfo'
    ]
)
