..
   Copyright (c) 2012-2015, Jeffrey N. Johnson
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

To date, Polymec has been a small project, with most of its present 
development done by a single person. This guide is not intended to be a 
comprehensive developer guide--it is only to establish best practices for 
maintaining the repository.

Repository Location
===================

The primary development repository for Polymec is presently 

https://bitbucket.org/jjphatt/polymec

Git Workflow
============

In the interest of simplicity, we have adopted the workflow advocated for 
by Peter Hintjens of ZeroMQ fame in http://hintjens.com/blog:24. This means 
that the main repository will only contain a master branch and that all 
development will occur in forks of this repository via pull requests. Forks 
may govern development as their owners see fit.

Disallowed File Types
=====================

In general, binary files should not be stored within the repository. These 
include (but are not limited to):

* Compiled assets of any sort.
* Documents (pdf), either from Polymec or in third-party libraries
* Binary data files (hdf5, silo, gz) within third-party libraries

