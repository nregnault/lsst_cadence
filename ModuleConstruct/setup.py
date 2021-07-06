from setuptools import setup, find_packages

import ALF_cadfunc

setup(
    name='ALF_cadfunc',
    version=ALF_cadfunc.__version__,
    packages=find_packages(),
    author='Angelo LF',
    description='Function for cadence metric',
    include_package_data=True,
    classifiers=[
        "Operati,g Systeme :: OS Independent",
        "Programming Language :: Python :: 3",
        ],
    license="WTFPL",
)
