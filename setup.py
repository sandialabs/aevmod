from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.test import test as TestCommand
from distutils.sysconfig import customize_compiler
import sysconfig
import glob
import os
import sys
import setuptools
import pybind11 as pb

base_path = os.path.dirname(__file__)

ext_modules = [
    # If you need to link extra libraries or specify extra include directories
    # see https://docs.python.org/3/extending/building.html#building-c-and-c-extensions-with-distutils
    Extension(
        'aevmod',
        glob.glob(os.path.join(base_path, 'src', '*.cpp')),
        include_dirs=[os.path.join(base_path,'include'),os.path.join(base_path,'include/sac'),os.path.join(base_path,'include/sac/mpl'),os.path.join(base_path,'include/sac/exptl'),pb.get_include()],
        language='c++',
        undef_macros=["NDEBUG"],
    ),
]

# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True

def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.

    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support is needed!')

class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        clang = False
        try:
            if os.environ['CC'] == "clang":
                clang = True
        except KeyError:
            pass
        if clang or sysconfig.get_config_vars()['CC'] == 'clang':
            c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):

        # remove this -Wstrict-prototypes flag for compiling c++ code
        customize_compiler(self.compiler)
        try:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
        except (AttributeError, ValueError):
            pass

        # add -std=c++11 or 14, and -fvisibility=hidden
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')

        for ext in self.extensions:
            ext.extra_compile_args = opts

        build_ext.build_extensions(self)

class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def run_tests(self):
        # import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)

setup(
    name='aevmod',
    version='1.0.0',
    author='Habib Najm',
    author_email='hnnajm@sandia.gov',
    long_description=open('README.rst').read(),
    ext_modules=ext_modules,
    cmdclass={
        'build_ext': BuildExt,
        'test': PyTest,
    },
    zip_safe=False,
    install_requires = [
        'numpy',
    ],
)
