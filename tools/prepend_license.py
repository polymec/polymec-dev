"""prepend_license.py -- prepends license information to all source files."""

license_text = """Copyright (c) 2012-2019, Jeffrey N. Johnson
All rights reserved.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/."""

old_license_text = """Copyright (c) 2012-2018, Jeffrey N. Johnson
All rights reserved.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/."""

# Here are the libraries we search for source files.
libs = ['core', 'geometry', 'solvers', 'model', 'io', 'mpi_serial']

import os, os.path

def find_sources(top_dir):
    sources = []
    for ldir in libs:
        libdir = os.path.join(top_dir, ldir)
        for root, dirs, files in os.walk(libdir):
            sources.extend([os.path.join(root, f) \
                for f in files if f.endswith('.h') or f.endswith('.h.in') or \
                                  f.endswith('.c') or f.endswith('.cpp') or \
                                  f.endswith('.f90.in')])
    return sources

def remove_old_license(source_file):
    f = open(source_file, 'r')
    lines = f.readlines()
    f.close()

    # Do we have the old license?
    old_license_lines = old_license_text.split('\n')
    has_old_license = True
    for i in range(len(lines)):
        if (i < len(old_license_lines)) and \
           (old_license_lines[i] not in lines[i]):
            has_old_license = False
            break

    f = open(source_file, 'w')
    if has_old_license:
        last_line = 1
        while old_license_lines[-1] not in lines[last_line]:
            last_line += 1
        for line in lines[last_line+1:]:
            f.write(line)
    else:
        for line in lines:
            f.write(line)
    f.close()

def add_new_license(source_file):
    f = open(source_file, 'r')
    lines = f.readlines()
    f.close()
    first_line = lines[0]
    f = open(source_file, 'w')
    license_lines = license_text.split('\n')
    for line in license_lines:
        if len(line) > 0:
            f.write('// ' + line + '\n')
        else:
            f.write('//\n')
    if license_lines[0] in first_line:
        last_line = 1
        while license_lines[-1] not in lines[last_line]:
            last_line += 1
        for line in lines[last_line+1:]:
            f.write(line)
    else:
        for line in lines:
            f.write(line)
    f.close()

import sys
if len(sys.argv) < 2:
    print('python prepend_license.py <dir>')
    exit(0)
root = sys.argv[1]

# Find the relevant source/header files.
sources = find_sources(root)

# Add the license.
for source in sources:
    if (old_license_text != ''):
        remove_old_license(source)
    add_new_license(source)
