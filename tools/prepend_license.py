"""prepend_license.py -- prepends license information to all source files."""

license_text = """Copyright (c) 2012-2015, Jeffrey N. Johnson
All rights reserved.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/."""

old_license_text = """Copyright (c) 2012-2015, Jeffrey N. Johnson
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this 
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."""

import os.path

def visit_files(sources, dirname, files):
    excluded_dirs = ['3rdparty', 'build']
    for ex_dir in excluded_dirs:
        if ex_dir in files:
            # Don't go there, girlfriend.
            files.remove(ex_dir)
    sources.extend([os.path.join(dirname, f) for f in files if f[-2:] == '.h'])
    sources.extend([os.path.join(dirname, f) for f in files if f[-2:] == '.c'])

def find_sources(dirname):
    sources = []
    os.path.walk(dirname, visit_files, sources)
    return sources

def remove_old_license(source_file):
    f = open(source_file, 'r')
    lines = f.readlines()
    f.close()
    first_line = lines[0]
    f = open(source_file, 'w')
    old_license_lines = old_license_text.split('\n')
    if old_license_lines[0] in first_line:
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
