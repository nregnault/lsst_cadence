# -*- mode: python; -*-

APPNAME = 'lsst_gaia_ubercal'
VERSION = '0.1.0'
top = '.'
out = 'build'
description = "Fast ubercal simulation software"

def options(opts):
    opts.load('python')

def configure(conf):
    conf.load('python')
    conf.check_python_version()

def build(bld):
    bld.recurse(['tools', 'py'])
