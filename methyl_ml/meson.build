project(
    'feature-vector-clustering',
    'cpp',
    version: '1',
    license: 'GPL-3',
    default_options: [
        'cpp_std=c++17',
        'warning_level=3',
        'optimization=3',
        'b_lto=true',
    ],
)

compiler = meson.get_compiler('cpp')

include_dirs = [
    include_directories('vendor/argparse/include'),
    include_directories('vendor/glob/single_include'),
    include_directories('vendor/indicators/include'),
    include_directories('vendor/libnpy/include'),
    include_directories('vendor/seqtk'),
]

dependencies = [
    dependency('openmp'),
    compiler.find_library('btllib'),
]

executable(
    'onehot',
    ['src/onehot.cpp'],
    dependencies: dependencies,
    include_directories: include_dirs,
    install: true,
)

executable(
    'bulk',
    ['src/bulk.cpp'],
    dependencies: dependencies,
    include_directories: include_dirs,
    install: true,
)