..
   Copyright (c) 2012-2018, Jeffrey N. Johnson
   All rights reserved.
   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

============================
Polymec Git Repository Guide
============================

This document describes the rules and workflows for working with Polymec's
Git repository.

Overview
========

So far, Polymec is a small project, with most development done by a single 
person. This guide isn't intended to be a comprehensive developer guide--it's 
mainly to establish best practices for maintaining the repository.

Repository Location
===================

The primary development repository for Polymec is 

https://github.com/polymec-dev/polymec-dev

Git Workflow
============

Polymec uses a forking workflow (https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow),
in which a contributor creates a fork of the primary repository, makes changes, 
and creates a pull request to start the merge process. The forking workflow 
is a simple and conservative flow that ensures that no changes are made to the 
primary repo without pull requests, during which code reviews, testing, and 
other diagnostics can be performed.

Disallowed File Types
=====================

In general, binary files aren't stored within the repository, since Git doesn't 
have a good strategy for tracking changes in them. Binary files include

* compiled assets of any sort.
* documents (pdf), either from Polymec or in third-party libraries
* binary data files (hdf5, silo, gz) within third-party libraries

Binary files that don't change may be added if doing so simplifies development.
