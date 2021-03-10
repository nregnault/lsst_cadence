from setuptools  import setup, find_packages

setup(
    name='sncadence', 
    version='0.1.0',
    author='Nicolas Regnault',
    author_email='nicolas.regnault@lpnhe.in2p3.fr',
    description='analysis of the sim/altsched/FBS cadences',
    include_package_data=True,
    packages=find_packages(),
    package_data={'html_report': ['data/html_report/*.js' 'data/html_report/*.css']}
)
