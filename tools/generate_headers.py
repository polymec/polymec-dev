# This script generates polymec.h and friends at the top level of the 
# build directory.
import sys, os.path
source_dir = sys.argv[1]
target_dir = sys.argv[2]

def visit_files(headers, dirname, files):
    if 'tests' in files:
        files.remove('tests')
    headers.extend([f for f in files if f[-2:] == '.h'])

def find_headers(dirname):
    headers = []
    os.path.walk(dirname, visit_files, headers)
    return headers

polymec_libs = ['core', 'geometry', 'integrators', 'model']

# Generate library-specific headers.
for lib in polymec_libs:
    headers = find_headers('%s/%s'%(source_dir, lib))
    header_name = '%s/polymec_%s.h'%(target_dir, lib)
    header = open(header_name, 'w')
    header.write('// polymec_%s.h -- automatically generated.\n'%lib)
    header.write('// This file is part of the polymec HPC library. See the license\n')
    header.write('// in the actual source files for details of distribution.\n\n')
    header.write('#ifdef __cplusplus\n')
    header.write('extern "C" {\n')
    header.write('#endif\n\n')
    if lib == 'core':
        header.write('#include "core/polymec.h"\n')
    for h in headers:
        header.write('#include "%s/%s"\n'%(lib, h))
    header.write('\n')
    header.write('#ifdef __cplusplus\n')
    header.write('}\n')
    header.write('#endif\n')
    header.close()

# Now generate the big guy.
header_name = '%s/polymec.h'%target_dir
header = open(header_name, 'w')
header.write('// polymec.h -- automatically generated.\n')
header.write('// This file is part of the polymec HPC library. See the license\n')
header.write('// in the actual source files for details of distribution.\n\n')
for lib in polymec_libs:
    header.write('#include "polymec_%s.h"\n'%lib)
header.write('\n')
header.close()

