"""prepend_license.py -- prepends license information to all source files."""

license_text = """Copyright 2012-2013 Jeffrey Johnson.

This file is part of Polymec, and is licensed under the Apache License, 
Version 2.0 (the "License"); you may not use this file except in 
compliance with the License. You may may find the text of the license in 
the LICENSE file at the top-level source directory, or obtain a copy of 
it at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License."""

old_license_text = ""

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
        while old_license_lines[-1] not in lines[i]:
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
    f.write('\n')
    if license_lines[0] in first_line:
        last_line = 1
        while license_lines[-1] not in lines[i]:
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
