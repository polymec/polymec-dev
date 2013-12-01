# This script generates polymec_version.h in the build directory.
import sys, subprocess, re
version_number = sys.argv[1]
header_file = sys.argv[2]

# Escape characters for use in C strings.
def escape(string, special_chars):
    for c in special_chars:
        if c == '"':
            string = string.replace('"', 'XXXDOUBLE_QUOTESXXX')
        elif c == '\n':
            string = string.replace('\n', '\\n"\n"')
        else:
            string = string.replace(c, '\\' + c)
    string = string.replace('XXXDOUBLE_QUOTESXXX', '\\"')
    return string

# Git revision ID
try:
    git_revision = ' (git revision %s'%subprocess.check_output(['git', 'log', '-1', '--format=format:%h']).strip()
    git_diff = subprocess.check_output(['git', 'diff']).strip()
    git_diff = escape(git_diff, ['"', '\n']) 
    if git_diff != '':
        git_revision += ', modified'
    git_revision += ')'
except:
    git_revision = ''
contents = """\n
// This file is automagically generated by update_polymec_version.py.
#ifndef POLYMEC_VERSION_H
#define POLYMEC_VERSION_H

static const char* POLYMEC_VERSION = "%s%s";
static const char* POLYMEC_GIT_DIFF = "%s";

#endif
"""%(version_number, git_revision, git_diff)
f = open(header_file, 'w')
f.write(contents)
f.close()

