"""prepend_license.py -- prepends license information to all source files."""

license_text = """Copyright (c) 2012-2017, Jeffrey N. Johnson
All rights reserved.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/."""

old_license_text = """Copyright (c) 2012-2016, Jeffrey N. Johnson
All rights reserved.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/."""

import os.path

def visit_files(sources, dirname, files):
    excluded_dirs = ['3rdparty', 'build']
    for ex_dir in excluded_dirs:
        if ex_dir in files:
            # Don't go there, girlfriend.
            files.remove(ex_dir)
    sources.extend([os.path.join(dirname, f) for f in files if f[-2:] == '.h'])
    sources.extend([os.path.join(dirname, f) for f in files if (f[-2:] == '.c' or f[-4:] == '.cpp')])

def find_sources(dirname):
    sources = []
    os.path.walk(dirname, visit_files, sources)
    return sources

def remove_old_license(source_file):
    f = open(source_file, 'r')
    lines = f.readlines()
    f.close()

    # Do we have the old license?
    old_license_lines = old_license_text.split('\n')
    has_old_license = True
    for i in xrange(len(lines)):
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
        f.write('// ' + line + '\n')
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

# Find the relevant source/header files.
sources = find_sources('.')

# Add the license.
for source in sources:
    if (old_license_text != ''):
        remove_old_license(source)   
    add_new_license(source)   
