#!/usr/local/bin/python

execfile("/usr/share/Modules/init/python.py")
module("purge")

# equivalent to source ~/.astrometry.src
module("load","gcc/4.9.1")
module("load","python/2.7.9")
module("load","cfitsio/3.370")
module("load","astrometry/0.50")

# equivalent to source ~/.sextractor.src
module("load","openmpi/1.6.1-gcc-4.4.6")
module("load","fftw/3.3.3-gcc-4.4.6-openmpi-1.6.1")
module("load","sextractor/2.8.6")
