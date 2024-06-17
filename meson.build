project('MethylScan', 'cpp', 'c',
        version : '0.0.1',
        license : 'GPL-3',
        default_options : [ 'cpp_std=c++17', 'warning_level=3', 'werror=true' ])


openmp_dep = dependency('openmp')

compiler = meson.get_compiler('cpp')
btllib_dep = compiler.find_library('btllib')
boost_math_dep = compiler.find_library('boost_math_c99')


add_project_arguments(
'-Wno-error=implicit-fallthrough=',
'-Wno-error=vla',
language : 'cpp'
)

add_project_arguments(
'-Wno-error=implicit-function-declaration',
language : 'c'
)



subdir('methyl_clus')
subdir('methyl_class')


rootpath = meson.source_root()
