from setuptools import setup
def readme():
    with open('README.rst') as f:
        return f.read()

setup(
    name = "niptplus",
    version = "0.1",
    description = "NIPTplus package",
    long_description = readme(),
    author = "bixichao",
    packages = ['niptplus'],
    include_package_data = True,
    zip_safe = False,
    test_suite = 'nose.collector',
    tests_require = ['nose'],
    scripts=['bin/nipt-joke'],
    entry_points = {
        'console_scripts': ['niptplus-commands=niptplus.command_line:niptplus_commands'],
    },
    install_requires = [
        'numpy',
        'scipy',
    ],
);

