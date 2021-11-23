#!/usr/bin/python3
from distutils.core import setup, Extension

module1 = Extension('const_l_gen',
                    sources = ['const_l_gen.cc'], 
		    extra_compile_args=['-std=c++11'])


setup (name = 'const_l_gen',
       version = '1.0',
       description = 'generate bases for n, l pairs',
       ext_modules = [module1])
