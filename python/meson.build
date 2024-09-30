project('yaeos', 'c',
  version : '0.1.0',
  license: 'MPL',
  meson_version: '>=0.64.0',
  default_options : ['warning_level=2'],
)

add_languages('fortran', native: false)

meson_path = meson.current_build_dir()

omp = dependency('openmp')
# lapack = dependency('lapack')

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()


incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()


incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()


# fpm
fflags = '-g -fPIC -funroll-loops -fstack-arrays -Ofast -frepack-arrays -faggressive-function-elimination -fopenmp'
cflags = '-fPIC'
fpm_build_dir = meson_path / 'fpm_build'

fpm = run_command(
    'fpm',
    'install',
    '--profile',
    'release',
    '--flag',
    fflags,
    '--c-flag',
    cflags,
    '--prefix',
    fpm_build_dir,
    check : false
)

# Error printing
message('fpm build output:')
output = fpm.stdout().strip()
message(output)
message('fpm build error:')
err = fpm.stderr().strip()
message(err)

# f2py
source_path = meson.global_source_root() / '..' / 'c_interface' / 'yaeos_c.f90'
c_wrapper = files(source_path)
module_name = 'yaeos_python'

f2py_gen_source = custom_target(module_name + '.c',
  input: c_wrapper,
  output: [module_name + 'module.c', module_name + '-f2pywrappers2.f90'],
  command: [
    py,
    '-m',
    'numpy.f2py',
    '@INPUT@',
    '-m',
    module_name,
    '--lower'
  ]
)

# =============================================================================
# Install Python sources
# -----------------------------------------------------------------------------
# yaeos root
python_lib = 'yaeos'
python_sources = [
    python_lib / '__init__.py',
    python_lib / 'core.py',
    python_lib / 'constants.py',
]
py.install_sources(python_sources, subdir:python_lib)


## lib
fapi_sources = [python_lib / 'lib' / '__init__.py']
py.install_sources(fapi_sources, subdir:python_lib / 'lib')


## models
models_sources = [python_lib / 'models' / '__init__.py']
py.install_sources(models_sources, subdir:python_lib / 'models')

### models.excess_gibbs
excess_gibbs_sources = [
  python_lib / 'models' / 'excess_gibbs' / '__init__.py',
  python_lib / 'models' / 'excess_gibbs' / 'nrtl.py',
  python_lib / 'models' / 'excess_gibbs' / 'unifac.py',
  python_lib / 'models' / 'excess_gibbs' / 'uniquac.py',
]
py.install_sources(excess_gibbs_sources, subdir:python_lib / 'models' / 'excess_gibbs')


### models.residual_helmholtz
residual_helmholtz_sources = [
  python_lib / 'models' / 'residual_helmholtz' / '__init__.py',
]
py.install_sources(residual_helmholtz_sources, subdir:python_lib / 'models' / 'residual_helmholtz')

#### models.residual_helmholtz.cubic_eos
cubic_eos_sources = [
  python_lib / 'models' / 'residual_helmholtz' / 'cubic_eos' / '__init__.py',
  python_lib / 'models' / 'residual_helmholtz' / 'cubic_eos' / 'cubic_eos.py',
  python_lib / 'models' / 'residual_helmholtz' / 'cubic_eos' / 'mixing_rules.py'
]
py.install_sources(cubic_eos_sources, subdir:python_lib / 'models' / 'residual_helmholtz' / 'cubic_eos')


# =============================================================================
# Extension fortran module
# -----------------------------------------------------------------------------
inc_np = include_directories(incdir_numpy, incdir_f2py, 'fpm_build' / 'include')

py.extension_module(module_name,
  [c_wrapper, f2py_gen_source],
  incdir_f2py / 'fortranobject.c',
  include_directories: [inc_np],
  dependencies : py_dep,
  link_args: ['-L' + 'fpm_build' / 'lib', '-l' + 'yaeos', '-llapack'],
  subdir: python_lib / 'lib',
  install : true,
)