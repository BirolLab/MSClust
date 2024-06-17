src_files = [ 'methyl_clus.cpp', 'city.cc']
inc_dirs = [ 'argparse' ]
deps = [ openmp_dep, btllib_dep ]
link_args = [ '-lbtllib' ]


executable('methyl-clus',
    src_files,
    dependencies : deps,
    include_directories : inc_dirs,
    install : true
)
