src_files = [ 'methyl_class.cpp', 'city.cc']
inc_dirs = [ 'argparse' ]
deps = [ openmp_dep, btllib_dep, boost_math_dep ]
link_args = [ '-lbtllib' ]


executable('methyl-class',
    src_files,
    dependencies : deps,
    include_directories : inc_dirs,
    install : true
)
