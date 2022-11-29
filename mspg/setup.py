import os
import sys
import sysconfig


setup_type = sys.argv[1]

#import setup and Extension functions from setuptools or distutils
#get numpy library directories

try:
    from setuptools import setup, Extension
    from setuptools.command.build_ext import build_ext
    use_setuptools = True
    print("setuptools is used.")

    # Trick to pip install numpy when it's not installed.
    # Reference: https://stackoverflow.com/questions/2379898/
    class CustomBuildExtCommand(build_ext):
        def run(self):
            import numpy
            self.include_dirs.append(numpy.get_include())
            build_ext.run(self)

except ImportError:
    from distutils.core import setup, Extension
    use_setuptools = False
    print("distutils is used.")

    try:
        from numpy.distutils.misc_util import get_numpy_include_dirs
    except ImportError:
        print("numpy.distutils.misc_util cannot be imported. Please install "
              "numpy first before installing spglib...")
        sys.exit(1)

# setup compile 
# Workaround Python issue 21121
config_var = sysconfig.get_config_var("CFLAGS")
if (config_var is not None and
    "-Werror=declaration-after-statement" in config_var):
    os.environ['CFLAGS'] = config_var.replace(
        "-Werror=declaration-after-statement", "")

sources = ['magnetic_cell.c',
           'magnetic_determination.c',
           'magnetic_primitive.c',
           'mspacegroup.c',     
           'magnetic_symmetry.c',
           'mspglib.c'
           ]


sources_spg = ['arithmetic.c',
               'cell.c',
               'delaunay.c',
               'debug.c',
               'determination.c',
               'hall_symbol.c',
               'mathfunc.c',
               'niggli.c',
               'overlap.c',
               'pointgroup.c',
               'primitive.c',
               'refinement.c',
               'sitesym_database.c',
               'site_symmetry.c',
               'spacegroup.c',
               'spg_database.c',
               'symmetry.c']

source_dir = './src'
source_dir_spg = './src/spglib'
include_dirs = [source_dir,source_dir_spg,]

if not use_setuptools:
    include_dirs += get_numpy_include_dirs()

for i, s in enumerate(sources):
    sources[i] = os.path.join(source_dir, s)

for i, s in enumerate(sources_spg):
    sources_spg[i] = os.path.join(source_dir_spg, s)

sources = sources + sources_spg

extra_compile_args = []
if setup_type == 'test':
    extra_compile_args.append("-UNDEBUG")
extra_link_args = []
define_macros = []

## Uncomment to activate OpenMP support for gcc
# extra_compile_args += ['-fopenmp']
# extra_link_args += ['-lgomp']

## For debugging
# define_macros = [('SPGWARNING', None),
#                  ('SPGDEBUG', None)]

# initialize Extension object
extension = Extension('mspglib._mspglib',
                      include_dirs=include_dirs,
                      sources=['_mspglib.c'] + sources,
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args,
                      define_macros=define_macros)


# Call setup function to build c extensions
if use_setuptools:
    setup(name='mspglib',
          version='1.0',
          cmdclass={'build_ext': CustomBuildExtCommand},
          setup_requires=['numpy', 'setuptools>=18.0'],
          description='This is a module to get magnetic spacegroup informations based on spglib.',
          author='suyl',
          install_requires=['numpy', ],
          packages = ['mspglib'],
          provides=['mspglib'],
          platforms=['all'],
          include_package_data=True,
          ext_modules=[extension])
else:
    setup(name='mspglib',
          version='1.0',
          description='This is a module to get magnetic spacegroup informations based on spglib.',
          author='suyl',
          requires=['numpy'],
          packages = ['mspglib'],
          provides=['mspglib'],
          platforms=['all'],
          include_package_data=True,
          ext_modules=[extension])
