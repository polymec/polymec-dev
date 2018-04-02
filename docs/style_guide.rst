..
   Copyright (c) 2012-2015, Jeffrey N. Johnson
   All rights reserved.
   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

===============================
Polymec Development Style Guide
===============================

This document describes rules, guidelines, and best practices for writing 
code for the Polymec library. The structure of this guide and parts of its 
prose were inspired by Google's C++ style guide, which appears at 

https://google-styleguide.googlecode.com/svn/trunk/cppguide.html

.. contents:: Table of Contents

Overview
========

Polymec is written in C. This language has come a long way since its inception, 
and its maintainers have given considerable effort to making C a desirable 
language for high-performance computing (HPC). Aside from the fact that all of 
the HPC libraries that Polymec uses are written in C (as most modern high-quality 
HPC scientific libraries have come to be), we find the following features of 
the language especially appealing:

* C is a simple, relatively small language, and it is usually clear to its 
  practitioners whether they understand what they are doing (or whether they 
  don't). Put another way, one can only get into so much trouble without asking 
  for it using this language.

* In C, what you get is what you ask for. There are very few hidden costs 
  that can be incurred accidentally.

* Since 1999, C has had many scientific programming features previously 
  available only in Fortran, such as dynamically-allocated multidimensional 
  arrays and complex numbers. Neither of these features is available in the 
  C++ core language.

* C is the common "hub" for several high-level and low-level programming 
  languages. Python, Ruby, C++, and Fortran all have standard mechanisms 
  for interfacing with it, for example.

* Link compatibility between C compilers removes much of the need for an 
  application to know "how polymec was built" (aside from specific 
  optimization settings).

* C still enjoys better support than C++ or Fortran within debuggers in UNIX 
  environments.

* Compile times for C are much shorter than for C++ or Fortran, often by 
  orders of magnitude.

In general, one can use any part of the C language described in the latest 
ISO standard, which at the time of this writing is the "C11" standard 
(ISO/IEC 9899:2011).

Source Code Organization
========================

Polymec comprises five C libraries:

* The ``core`` library, which contains utilities, scientific data structures,
  and components that can be used in applications. This library is used by the 
  other three Polymec libraries.

* The ``geometry`` library, which consists of various implementations of 
  algorithms for mesh generation, implicit functions and space-time functions, 
  and geometric operations in general. This library uses the ``core`` library.

* The ``solvers`` library, which has classes and tools for solving 
  nonlinear equations. This library uses the ``core`` library.

* The ``model`` library, which contains high-level data structures for 
  constructing science applications and numerical models. This library uses 
  the ``core`` and ``geometry`` libraries.

* The ``io`` library, which contains high-performance input/output classes
  and functions that support most of the important types from the other 
  libraries. This library uses the ``core``, ``model``, and ``geometry`` 
  libraries.

Each Polymec source file belongs to one of these libraries, and is located 
in the subdirectory of Polymec named after that library.

Header Files
============

In general, there should be a header file for each signicant type or "class" 
in Polymec. A header file for a class named a_class would be named a_class.h.
In some cases, a single function (unrelated to a type) may occupy a header 
file, and that header file would be named after the function.  In others, a 
header file may contain a set of related functions, and its name should 
concisely reflect the purpose of those functions.

Self-contained Headers 
----------------------

Header files should be self-contained and have a .h suffix. In other words, 
it should be possible to include a header file without regard for rules 
relating to the order of its inclusion, or for other headers that are 
"understood" to be included when it is used. Specifically, a header file 
should use header guards, should include all the files that it needs, and 
should not require any particular symbols to be defined.

Header Guards
-------------

A header file should have #define guards to prevent multiple inclusion. The 
format of the guard should be ``POLYMEC_<HEADER_NAME>_H``, e.g.:

``#ifndef POLYMEC_FOO_H
#define POLYMEC_FOO_H``

"C++" guards that use the extern "C" specification are not necessary, since 
Polymec generates higher-level headers that are safe for inclusion in 
C++ programs.

Including Headers Within a Polymec Header 
-----------------------------------------

Any other header files included in the header should be included in the 
following order:

1. System-level headers
2. Third-party library headers
3. Polymec library headers.

System-level headers must be enclosed in angle brackets, while third-party 
library headers and Polymec headers should be enclosed in double quotes.

When including Polymec header files in a library that belongs to the Polymec
library itself, one must specify the full path of the header file relative to 
the top level Polymec source directory in the ``#include`` directive, 
e.g.

``#include "core/polymec.h``

Classes 
-------

A "class" in Polymec is a struct with a set of associated functions. The 
struct must only be declared in its header file--its body should NOT be 
defined in a header file unless doing so is required for technical reasons. 
The body should be defined in the source (.c) file associated with the header.
Defining class bodies in source files is a common practice in C to reduce 
code coupling, and it is emulated in the C++ "Pimpl" idiom.

Functions 
---------

Any function that is part of Polymec's API should be declared within a 
header file. A function may be "inlined" using the ``static inline``
C construct. Functions with no arguments must be declared with ``void`` in 
their argument list, per the C11 standard.

Global variables 
----------------

No global variables should appear within a header file, apart from constants 
(which are preferred to macros, since they can be checked by the compiler). 
Mutable global variables should be restricted to translation units in which 
they are manipulated. If a global data structure needs to exist, an appropriate 
interface should be designed and implemented in terms of functions.

Other Symbols 
-------------

Inlined functions should be used instead of macros where possible. Similarly, 
constants should be used instead of macros where possible.

Special Datatypes
=================

In Polymec, floating point variables should be stored using the 
``real_t`` type. Integers representing indices that can assume 
large values should be stored using the ``index_t`` type. Both of 
these types are declared in ``core/polymec.h``.

Polymec's ``core`` library contains several standard data structures such as 
dynamic arrays, tuples, linked lists, unordered and ordered sets, unordered 
and ordered maps, tables, space-filling curves, space-time functions, sparse 
matrices, random number generators, kd-trees, and so on. Please check here 
before you decide to implement your own data structure.

Floating point comparisons
--------------------------

Avoid comparing two floating point numbers to see whether they are equal. 
Instead, use the ``reals_equal`` and ``reals_nearly_equal`` functions that 
allow two floating point numbers to be equated when they differ by an amount 
below some threshold. The former uses a threshold set by ``set_real_epsilon``, 
whereas the latter uses a threshold specified in the call.

Scoping
=======

Static functions 
----------------

A function that is used only within one translation unit should be declared 
static so that its name does not appear in the list of symbols for the 
Polymec library.

Local variables 
---------------

A local variable should be declared as close as possible to the location(s) 
at which it is used. This makes it easier to identify problems involving 
that variable.

A variable should be initialized where it is declared, unless such an 
initialization renders a code construction awkward or inefficient.

Scoping operators
-----------------

If a function has a large number of localized variables that perform work, 
curly braces should be used to create a local scope containing these variables.
This eases the process of debugging functions by eliminating these variables 
from portions of the function that don't use them.

Classes
=======

As mentioned in the section on header files, a Polymec class consists of a 
struct representing that class, and an associated set of functions that
are considered its methods. Class bodies are defined in source files 
only, unless their internal structure is intended to be explicitly exposed to 
developers. A class type should be "typedefed" so the 
``struct`` keyword is not required to precede it.

The struct and functions defining a class are governed by the following set of 
conventions.

Class type (struct) 
-------------------

The struct representing the class type should end in ``_t``. For 
example, if we declare a "point" class, we might declare a struct

``typedef struct point_t;``

in a header file (point.h, say), and define the struct in a source file 
(e.g. point.c).

Class constructor(s)
--------------------

Typically, a class will have a single constructor function named 
``<CLASS>_new`` that takes a number of arguments for initializing the class, 
and returns a newly-allocated pointer to an instance of the corresponding 
class struct. For example, we might define a constructor for our point class 
thus:

``point_t* point_new(real_t x, real_t y, real_t z);``

Sometimes more than one constructor will be necessary, or a constructor that 
converts another datatype to a given instance of a class will be convenient.
In this case, each constructor should briefly convey its nature. For example, 
a constructor that converts an array of ``real_t`` to a point might 
be declared 

``point_t* point_from_array(real_t* array);``

A constructor function should take any arguments it needs to completely 
initialize an variable of that class type, and return a pointer to such an 
initialized variable. We refer to these variables as objects.

Class destructor 
----------------

A single destructor function must be defined for any class that does not use 
garbage collection. The destructor function must have no return type, and must 
be named ``<CLASS>_free``. It must take as an argument a pointer 
to the struct representing an instance of that class. For example:

``void point_free(point_t* point);``

The destructor must completely deallocate any resources allocated to the 
argument object in its construction and during its lifetime.

Methods 
-------

A method for a class should be named ``<CLASS>_<METHOD>`` and should 
always take a pointer to the struct representing an instance of 
that class as its first argument. For example, the following method returns 
the distance between the given point and another point:

``real_t point_distance(point_t* point, point_t* other);``

Methods should be defined in a manner similar to the idioms found in modern 
object-oriented programming languages such as C++ and Java. After the first 
argument, arguments should be ordered with input values at the beginning 
of the argument list followed by output values at the end.

Polymorphism in C 
-----------------

Polymorphic classes in Polymec have "abstract" base classes with virtual 
tables that dispatch calls to functions in the class interface. The base class 
consists of:

1. A class type struct possessing a context pointer for an instance
2. A virtual table (vtable) struct consisting of a set of function pointers 
   matching the interface for the class
3. A constructor function that creates a descendant object using a context 
   pointer, a vtable, and any other data needed.
4. Any other functions needed to implement a destructor and/or methods for the 
   polymorphic class.

This approach to polymorphism is called "prototype polymorphism," and is used 
in some other programming languages such as Lua. The idea is that the behavior 
of a polymorphic class is tied to a specific instance of that class, not to its 
type. 

One virtue of this approach is that a single "object" (represented by a 
context pointer) can assume many different roles as a subtype of several 
base classes, using several different virtual tables. In a sense, this 
ability resembles that of the ``interface`` idiom in the Java and C# 
programming languages, avoiding the difficulties of multiple inheritance one 
encounters in C++.

See Polymec's ``model`` class in ``model/model.h`` and ``model/model.c`` for 
an example of how polymorphic data structures can be implemented using this 
model.

Structs as "Plain Old Datatypes" (PODs)
=======================================

Occasionally, it may be expedient to declare a struct representing a simple 
container, or "Plain Old Datatype" (POD). In this case, no constructor or 
destructor or methods are needed for manipulation unless such mechanisms make 
the POD more convenient to use.

Functions
=========

Functions not associated with classes follow very similar guidelines to 
methods: input arguments come before output arguments.

Length of a Function Body
-------------------------

There is no formal limit to the length of a Polymec function implementation. 
If breaking up a function into separate functions is practical, feel free to 
do so. However, creating lots of ancillary structure just to break up a long 
function is counterproductive. Use your judgement.

The function indeed may be poorly designed if it is difficult to break up. 
On the other hand, if the function is performing a complicated task with lots 
of tightly coupled steps, attempting to break it up may further obfuscate its 
task.

At the end of the day, arguments about the optimal length of a function are 
based in aesthetics and often exert strange and unnatural pressures on 
code development, encouraging people to write code with few comments, lots of 
side effects, and/or excessive numbers of tightly-coupled "sub-functions."

Memory Management
=================

Memory is typically allocated in polymec by calling the ``polymec_malloc`` 
function. This function has the same signature as the standard C ``malloc`` 
function, but allows access to some of Polymec's special memory allocation 
features.

Ownership
---------

To minimize complexity, try to assign a single owner to an allocated resource. 
Try to avoid ownership transfers, as these can create complicated resource 
management issues. In typical HPC programming patterns, ownership transfers 
are not usually necessary for objects using large amounts of resources.

When an object needs to store dynamically-allocated data (e.g. from an array 
or another object) that is assigned to it via a method, that data should 
typically be copied to an array within the object. An alternative to copying 
the data is to "transfer ownership" to the object, in which case the object 
becomes responsible for deallocating the data when it is destroyed. In this 
case, it is sometimes said that the object "consumes" the data, meaning that 
the data cannot be assumed to be available after the object is deleted, and 
that the integrity of the data cannot be guaranteed after the transfer of 
ownership.

Polymec does not use smart pointers. Smart pointers require special machinery 
to achieve thread safety, and often do not perform well when the number of 
threads in a process becomes significant. Moreover, the excessive use of 
smart pointers can lead to complexities in the relationships of objects to 
one another, making a system difficult to analyze.

If two or more objects need access to a resource that is not clearly owned by 
one or the other, try to consider the relationship between these objects and 
determine whether the resource should be owned by one of these objects or by 
another "service" object. If the resource is small and/or it can safely be 
destroyed in a non-deterministic fashion after it has been used, consider 
managing it with garbage collection.

In any case, any functions/methods that perform ownership transfers should 
describe the transfer in their documentation.

Garbage Collection
------------------

Classes representing small objects whose ownership is not clear-cut may use 
a simple reference counting system supported by Lua. An object of a 
garbage-collected type has no destructor in its API, since its destruction is 
performed automatically some time after all references to it have been 
destroyed.

Use ``polymec_gc_malloc`` to allocate garbage-collected resources. Any 
resources managed by the garbage-collected object should be allocated with 
``polymec_malloc`` and freed with ``polymec_free``. The resources should be 
freed in a private destructor supplied to ``polymec_gc_malloc``. 

In Lua, garbage-collected objects are managed automatically. In C, you need 
to increment the reference count of an object you need to keep around. Do 
this with ``polymec_retain``. When you're finished with an object and you 
don't care when it gets collected, release it with ``polymec_release``.

For an example of a simple garbage-collected type in Polymec, see the ``point``
class in ``core/point.h`` (and its implementation in ``core/point.c``). For a 
more elaborate example of a garbage-collected type that manages its own 
resources, see the ``st_func`` class in ``core/st_func.h`` and 
``core/st_func.c``.

Special Allocators
------------------

Polymec has a few specialized allocators intended to reduce the overhead of 
requesting memory from the operating system. These are:

* ``std_allocator`` - The standard C allocator ``malloc``. This is used by default.

* ``arena_allocator`` - An allocator that preallocates a large "arena" from which 
  memory is dispensed to requestors. The arena has a large initial size and then 
  is resized as necessary. Like the operating system's heap, this arena can 
  become fragmented over time if memory allocations are unstructured. However,
  the number of memory requests to the operating system is minimized by 
  pre-allocating the arena beforehand.

* ``pool_allocator`` - An allocator that dispenses memory in several "pools". 
  Essentially, this is a segmented version of the ``arena_allocator``, and can 
  be more flexible.

Allocators are controlled via an interface defined in ``core/allocators.h``. 
You will need to experiment with these allocators to gain an understanding of 
their benefits, drawbacks, and general capabilities.

Naming
======

Names of structs, classes, and enumerated types should all contain only 
lower-case characters with words separated by underscores, ending in 
``_t``. Abbreviations are allowed if their meaning is reasonably clear. For 
example: ``mesh_t``, ``point_t``, ``ode_integrator_t``. ``adj_graph_t``.

Function and method names should also use only lower-case letters with 
words separated by underscores. Unintelligible abbreviations should not be 
used for struct, class, or function names. Examples of good function and 
method names are ``point_distance``, ``partition_mesh``, and 
``polymec_timer_start``.

Similarly, a variable (local or global, including fields in structs and classes)
should strive to use only lower-case letters with words separated by 
underscores. Exceptions can be made if it makes code clearer. For example, 
capital letters and/or abbreviations may help a variable representing a 
quantity resemble a mathematical symbol whose role is clear from the context 
in which it is used. Use your judgement. Examples of good variable names are 
``adj_graph``, ``mesh``, ``model``, ``precond``, ``integ``, and ``xc``.

Constants, fields within enumerated types, and preprocessor macros should use 
all capital letters with words separated by underscores. If these appear in 
header files, they should have descriptive names that are unique within the 
library. An example of a good name for a constant is ``REAL_EPSILON``. 
Examples of enumerated type fields are ``MESH_NODE``, ``MESH_EDGE``, 
``MESH_FACE``, and ``MESH_CELL``. Examples of preprocessor macros are 
``START_FUNCTION_TIMER``, ``DECLARE_2D_ARRAY``, and ``DEFINE_UNORDERED_SET``.

Comments
========

Use C++ style comments (``//``), which have been supported in C since the 
C99 standard. C-style comments (``/* */``) are clunkier and harder for 
editors to parse correctly.

Polymec does not use an inline documentation generator like Doxygen, so there 
are no documentation markup tags. Class types, structs, and enumerated types 
should be commented with a brief synopsis of their purpose. The comments 
should precede the ``typedef`` for the type.

A function should be commented with a brief description of the function, its 
preconditions, postconditions, and return values where applicable. The 
comments should precede the function declaration in header files.

It is not necessary to repeat comments for a class/function in its 
implementation file if its declaration hsa been documented. However, classes 
and functions that appear only within an implementation file should be 
documented for completeness.

Formatting
==========

The following formatting rules are non-negotiable for source code in Polymec:

* Use 2 spaces per indentation level.
* No tabs are allowed in source files -- use only spaces.

The following guidelines are offered for readably formatted code:

* If a function declaration doesn't fit neatly on a line, break the line after 
  an argument and align the following argument with its first. As long as the 
  declaration and definition are clearly readable, it's fine.
* Curly braces that open and close new scopes each go on their own line, not 
  at the end of a line containing other code.
* If a line is excessively long (in other words, if it doesn't fit on a single 
  screen on a luxuriously large monitor), consider breaking it up.
* C preprocessor directives are not indented at all.
* For functions with several parameters, consider linebreaks after each 
  parameter, and consider aligning the parameters to improve readability.
