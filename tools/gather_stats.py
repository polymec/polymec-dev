"""gather_stats.py -- gathers statistics for the polymec code base."""

import os.path

def visit_files(headers_and_sources, dirname, files):
    excluded_dirs = ['3rdparty', 'build']
    for ex_dir in excluded_dirs:
        if ex_dir in files:
            # Don't go there, girlfriend.
            files.remove(ex_dir)
    headers_and_sources[0].extend([os.path.join(dirname, f) for f in files if f[-2:] == '.h'])
    headers_and_sources[1].extend([os.path.join(dirname, f) for f in files if f[-2:] == '.c'])

def find_sources(dirname):
    headers_and_sources = ([], [])
    os.path.walk(dirname, visit_files, headers_and_sources)
    return headers_and_sources

def count_source_lines(source_files):
    num_lines = 0
    num_lines_wo_comments = 0
    for f in source_files:
        fp = open(f, 'r')
        lines = fp.readlines()
        for line in lines:
            line = line.strip()
            if line != '':
                num_lines += 1
                if not line.startswith('//'):
                    num_lines_wo_comments += 1
        fp.close()
    return (num_lines, num_lines_wo_comments)

def count_types_and_functions(headers):
    num_types = 0
    num_functions = 0
    for f in headers:
        in_func_decl = False
        fp = open(f, 'r')
        lines = fp.readlines()
        for line in lines:
            line = line.strip()
            if not line.startswith('//') and not line.startswith('#'):
                if 'typedef struct' in line:
                    num_types += 1
                elif 'typedef' not in line:
                    if '(' in line and ')' in line and ';' in line:
                        num_functions += 1
                    elif '(' in line and not ';' in line:
                        in_func_decl = True
                    elif not '(' in line and ')' in line and ';' in line and in_func_decl:
                        in_func_decl = False
                        num_functions += 1
        fp.close()
    return (num_types, num_functions)

def count_others(source_files):
    num_fixmes = 0
    for f in source_files:
        fp = open(f, 'r')
        lines = fp.readlines()
        for line in lines:
            line = line.strip()
            if 'FIXME' in line:
                num_fixmes += 1
        fp.close()
    return num_fixmes

# Find the relevant source/header files.
headers, sources = find_sources('.')
num_header_files = len(headers)
num_source_files = len(sources)

# Count the total number(s) of source lines.
(num_lines, num_lines_wo_comments) = count_source_lines(headers + sources)

# Count function points in headers.
(num_types, num_functions) = count_types_and_functions(headers)

# Miscellany.
num_fixmes = count_others(headers + sources)

print '-----------------------'
print 'Statistics for polymec:'
print '-----------------------'
print 'number of header files:                    %i'%num_header_files
print 'number of source files:                    %i'%num_source_files
print 'number of source lines (without comments): %i (%i) '%(num_lines, num_lines_wo_comments)
print 'number of types:                           %i'%num_types
print 'number of functions:                       %i'%num_functions
print 'number of FIXMEs:                          %i'%num_fixmes

