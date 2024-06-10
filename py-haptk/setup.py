
from setuptools import find_packages, setup

setup(
    name='haptk',
    # packages=find_packages(include=['haptk']),
    packages=find_packages(),
    version='0.0.1',
    description='HAPTK Python library',
    author="xosxos",
    install_requires=['json', 'gzip', 'ete3', 'rustworkx', 'argparse', 'pyqt5'],
    setup_requires=['pytest-runner'],
    tests_require=['pytest==4.4.1'],
    test_suite='tests',
)
