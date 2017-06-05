# This script generates polymec.h and friends at the top level of the 
# build directory.
import sys, os.path, shutil
source_dir = sys.argv[1]
target_dir = sys.argv[2]

def visit_c_files(headers, dirname, files):
    if 'tests' in files:
        files.remove('tests')
    h_files = [f for f in files if f[-2:] == '.h' and f[0] != '.']
    headers.extend(h_files)

def find_c_headers(dirname):
    headers = []
    os.path.walk(dirname, visit_c_files, headers)
    return headers

def files_match(file1, file2):
    try:
        f1 = open(file1, 'r')
        f2 = open(file2, 'r')
        result = f1.read() == f2.read()
        f1.close()
        f2.close()
    except:
        result = False
    return result

polymec_libs = ['core', 'geometry', 'solvers', 'model', 'meshless']

# Generate library-specific headers.
for lib in polymec_libs:
    c_headers = find_c_headers('%s/%s'%(source_dir, lib))
    header_name = '%s/polymec_%s.h'%(target_dir, lib)
    staging_header_name = header_name + '.new'
    header = open(staging_header_name, 'w')
    header.write('// polymec_%s.h -- automatically generated.\n'%lib)
    header.write('// This file is part of the polymec HPC library. See the license\n')
    header.write('// in the actual source files for details of distribution.\n\n')
    LIB = lib.upper()
    header.write('#ifndef POLYMEC_%s_LIBRARY_H\n'%LIB)
    header.write('#define POLYMEC_%s_LIBRARY_H\n\n'%LIB)
    header.write('#ifdef __cplusplus\n')
    header.write('extern "C" {\n')
    header.write('#endif\n\n')
    if lib == 'core':
        header.write('#include "core/polymec.h"\n')
    for h in c_headers:
        header.write('#include "%s/%s"\n'%(lib, h))
    header.write('\n')
    header.write('#ifdef __cplusplus\n')
    header.write('}\n')
    header.write('#endif\n\n')
    header.write('#endif\n\n')
    header.close()

    # Check to see whether the new header's contents match the old one's 
    # identically. If they don't match, copy the new one over.
    if not files_match(header_name, staging_header_name):
        shutil.copyfile(staging_header_name, header_name)
    os.remove(staging_header_name)

# Now generate the big guy.
header_name = '%s/polymec.h'%target_dir
staging_header_name = header_name + '.new'
header = open(staging_header_name, 'w')
header.write('// polymec.h -- automatically generated.\n')
header.write('// This file is part of the polymec HPC library. See the license\n')
header.write('// in the actual source files for details of distribution.\n\n')
header.write('#ifndef POLYMEC_LIBRARY_H\n')
header.write('#define POLYMEC_LIBRARY_H\n\n')
for lib in polymec_libs:
    header.write('#include "polymec_%s.h"\n'%lib)
header.write('\n')
header.write('#endif\n\n')
header.close()

# Check to see whether the new header's contents match the old one's 
# identically. If they don't match, copy the new one over.
if not files_match(header_name, staging_header_name):
    shutil.copyfile(staging_header_name, header_name)
os.remove(staging_header_name)
