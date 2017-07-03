[![Build Status](https://travis-ci.org/polymec/polymec-dev.svg?branch=master)](https://travis-ci.org/polymec/polymec-dev)
[![Coverage Status](https://coveralls.io/repos/github/polymec/polymec-dev/badge.svg?branch=master)](https://coveralls.io/github/polymec/polymec-dev?branch=master)

Polymec
=======

Copyright (c) 2012-2017, Jeffrey N. Johnson
All rights reserved.

Polymec is a set of libraries that can be used to construct models for 
physical systems. These libraries are designed to allow a computational 
scientist to construct low-maintenance science applications with robust 
capabilities. 

Features
--------

* A comprehensive set of C tools appropriate for high performance computing, 
  including containers, spatial data structures, and polished interfaces to 
  important libraries and facilities.
* A layered structure allowing different levels of buy-in, from using only 
  basic data structures to adopting a framework for building high performance 
  simulators.
  represented by unstructured (polyhedral) grids and/or mesh-free methods.
* Stiffly-accurate, implicit time integrators that use high-quality 
  preconditioners without a lot of application-specific classes.
* A strategy that favors a small code footprint and good scalability, likely 
  at the expense of implementations that run optimally on serial computers.

License
-------

Polymec is licensed under the Mozilla Public License (MPL) version 2.0, which 
is defined in the LICENSE file, and also available at http://mozilla.org/MPL/2.0/.
This license allows you to use Polymec in a commercial code as long as 
certain requirements are satisfied. Please read or read about the license.

The 3rd-party libraries used by Polymec have their own licenses which are 
described within the various source trees in the 3rdparty directory. Source 
files have been minimally modified within these trees; typically they contain 
files extracted directly from distribution tarballs. Instructions for 
obtaining the distributions are given within.

