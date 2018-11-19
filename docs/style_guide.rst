..
   Copyright (c) 2012-2018, Jeffrey N. Johnson
   All rights reserved.
   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

===============================
Polymec Development Style Guide
===============================

This document describes rules, guidelines, and best practices for writing 
code for the Polymec library. The structure of this guide was inspired by 
Google's C++ style guide, which appears at 

https://google.github.io/styleguide/cppguide.html

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
  languages. Python, Ruby, Lua, C++, and Fortran all have standard mechanisms 
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

Polymec consists of five C libraries:

* The ``core`` library, which contains utilities, scientific data structures,
  and components that can be used in applications. This library is used by the 
  other libraries.

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

Self-Contained Headers 
----------------------

Header files are self-contained and have a .h suffix. A "self-contained" header
file can be included in a translation unit without regard for rules 
relating to the order of its inclusion, or for other headers that are 
"understood" to be included when it is used. 

Briefly, a header file
* requires header guards
* should include all the files that it needs
* should not require any particular symbols to be defined

Header Guards
-------------

A header file uses #define guards to prevent multiple inclusion. The 
format of the guard is ``POLYMEC_<HEADER_NAME>_H``, e.g.:

``#ifndef POLYMEC_FOO_H
#define POLYMEC_FOO_H``

"C++" guards that use the extern "C" specification are not necessary, since 
Polymec generates higher-level headers that are safe for inclusion in 
C++ programs.

Including Headers within a Polymec Header 
-----------------------------------------

Any other header files included in the header should be included in the 
following order:

1. System-level headers
2. Third-party library headers
3. Polymec library headers.

System-level headers are enclosed in angle brackets, while third-party 
library headers and Polymec headers are enclosed in double quotes.

When including Polymec header files in a library that belongs to the Polymec
library itself, one specifies the full path of the header file relative to 
the top level Polymec source directory in the ``#include`` directive, 
e.g.

``#include "core/polymec.h``

Classes and Structs
-------------------

A "class" in Polymec is a struct with a set of associated functions. The 
struct must only be declared in its header file--its body should NOT be 
defined in a header file unless doing so is required for technical reasons. 
The body should be defined in the source (.c) file associated with the header.
Defining class bodies in source files is a common practice in C to reduce 
code coupling, and it is emulated in the C++ "Pimpl" idiom.

In the context of Polymec, a "struct" is a C struct containing data that can 
be freely exposed. Structs have no behavior and internal state to manage. They 
are defined in header files so that their data members are visible and accessible.
Examples of structs are the ``point_t`` and ``vector_t`` types, which represent
points and vectors in 3D space.

Functions 
---------

Any function that is part of Polymec's API is declared within a header file. 
A function may be "inlined" using the ``static inline`` C construct. Functions 
with no arguments are declared with ``void`` in their argument list, in 
accordance with the C11 standard.

Global variables 
----------------

Avoid global variables in header files, apart from constants (which are 
preferred to macros, since they can be checked by the compiler). Mutable 
global variables should be restricted to translation units in which they are 
manipulated, and should be declared as ``static``. If you must expose a global 
resource, design an appropriate interface so that it can be properly managed.

Other Symbols 
-------------

Use inlined functions instead of macros where possible. Similarly, use 
constants instead of macros where possible.

Special Datatypes
=================

Polymec does a lot of stuff with real numbers, and sometimes even with 
complex and imaginary ones. Use the ``real_t`` type to represent real 
quantities, and the ``complex_t`` type. You can pass real numbers around with 
MPI using ``MPI_REAL_T``, and complex numbers with ``MPI_COMPLEX_T``. These 
types are defined properly for the level of precision for which Polymec is 
configured.

You might also want a 64-bit index type if you're assembling a distributed 
linear system or a giant adjacency graph. The ``index_t`` type is guaranteed 
to be 64-bit, and you can use the MPI type ``MPI_INDEX_T`` to pass it around 
between processes.

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

Static Functions 
----------------

A function that is used only within a single translation unit should be 
declared static within that translation unit. This prevents its name from 
appearing in the list of exported symbols for the Polymec library.

Local Variables 
---------------

Declare a local variable as close as possible to where it is used. This makes 
it easier to identify issues involving that variable.

Initialize a variable when you declare it wherever practical.

Scoping Operators
-----------------

If a function has a large number of localized variables that perform work, 
curly braces should be used to create a local scope containing these variables.
This eases the process of debugging functions by eliminating these variables 
from portions of the function that don't use them.

Classes
=======

As mentioned in the section on header files, a Polymec class consists of a 
struct representing that class, and an associated set of functions that
are considered its methods. Define class bodies in source files only, unless 
their internal structure is intended to be explicitly exposed to developers. 
"Typedef" your class type so the ``struct`` keyword can be omitted from its 
type.

The struct and functions defining a class are governed by a few simple 
conventions.

Class Type (Struct) 
-------------------

The struct representing the class type ends in ``_t``. For example, if you 
declare a "washing machine" class, you might declare a struct

``typedef struct washing_machine_t;``

in a header file (``washing_machine.h``, say), and define the struct in a 
source file (e.g. ``washing_machine.c``).

Class Constructor(s)
--------------------

Typically, a class has a single constructor function named ``<CLASS>_new`` 
that takes a number of arguments for initializing the class, and returns a 
newly-allocated pointer to an instance of the corresponding class struct. 
For example, we might define a constructor for our point class thus:

``point_t* point_new(real_t x, real_t y, real_t z);``

Sometimes it's convenient to provide more than one constructor, or a 
constructor that converts another datatype to a given instance of a class. 
In these cases, name each constructor so that it briefly conveys its purpose. 
For example, a constructor that converts an array of ``real_t`` to a point 
might be declared 

``point_t* point_from_array(real_t* array);``

A constructor function takes any arguments it needs to completely initialize 
an variable of that class type, and returns a pointer to such an initialized 
variable. We refer to these variables as objects.

Class destructor 
----------------

Define a single destructor function for any class that does not use reference 
counting. The destructor function has no return type, and must be named 
``<CLASS>_free``. The destructor take a single argument: a pointer to the 
struct representing an instance of that class. For example:

``void point_free(point_t* point);``

The destructor completely deallocates any resources allocated to the 
argument object in its construction and during its lifetime.

Methods 
-------

Name a method for a class using that class's name as a prefix: 
``<CLASS>_<METHOD>``. The first argument to a method is a pointer to the 
struct representing an instance of that class. For example, the following 
method returns the distance between the given point and another point:

``real_t point_distance(point_t* point, point_t* other);``

Define methods just as you would in contemporary object-oriented programming 
languages like C++ and Java. If it's practical, lead the list of parameters 
with input values, and put output parameters at the end.

Polymorphism in C 
-----------------

Polymorphic classes in Polymec have "abstract" base classes with virtual 
tables that dispatch calls to functions in the class interface. The base class 
consists of:

1. A class type struct possessing a context pointer for an instance
2. A virtual table (vtable) struct consisting of a set of function pointers 
   matching the interface for the class
3. A constructor function that creates a descendant object using a context 
   pointer, a vtable, and any other data needed
4. Any other functions needed to implement a destructor and/or methods for the 
   polymorphic class

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

Sometimes it's convenient to declare a struct representing a simple container, 
or "Plain Old Datatype" (POD). In this case, no constructor or destructor or 
methods are needed for manipulation unless such mechanisms make the POD more 
convenient to use.

Functions
=========

Functions not associated with classes follow very similar guidelines to 
methods: input arguments come before output arguments.

Length of a Function Body
-------------------------

There is no formal limit to the length of a Polymec function implementation. 
If breaking up a function into separate functions is practical, feel free to 
do so. However, creating lots of ancillary structure just to break up a long 
function can be counterproductive. Use your judgement.

A function may be poorly designed if it is difficult to break up. On the other 
hand, if the function performs a complicated task with lots of tightly-coupled 
steps, attempting to break it up may make it even more confusing.

At the end of the day, arguments about the optimal length of a function are 
aesthetic. These arguments often exert strange and unnatural pressures on code 
development. At worst, they encourage people to write code with few comments, 
lots of side effects, and/or excessive numbers of tightly-coupled 
"sub-functions." Your mileage may vary.

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

Reference Counting
------------------

Classes representing small objects whose ownership is not clear-cut may use 
a simple reference counting system supported by Lua. An object of a 
reference-counted type has no destructor in its API, since its destruction is 
performed automatically some time after all references to it have been 
destroyed.

Use ``polymec_refcounted_malloc`` to allocate a reference-counted resource. 
Allocate any resources managed by this reference-counted object as usual, with 
``polymec_malloc``. Free these resources in with ``polymec_free`` in a private 
destructor you supply to ``polymec_refcounted_malloc``. 

A newly-created reference-counted object starts with a reference count of 1. 
If you want to retain a reference to an object ``o`` to prevent it from being 
destroyed, increment its reference count like this:

``retain_ref(o);``

If you don't need to bump the reference count, but you want to annotate your 
code to say you're "borrowing" the object ``o`` temporarily, you can be a good 
citizen and use 

``borrow_ref(o);``

This statement doesn't actually do anything--it's just good manners. Finally, 
when you're done with the reference-counted object, you can release it with 

``release_ref(o);``

Objects can be created in C or in the Lua interpreter, and these objects 
can be transferred from one environment to the other. Lua maintains its own 
resources with automatic garbage-collection and will respect the reference 
count of resources you transfer to it. See ``lua_transfer_object`` for details
on how to migrate an object one way or the other.

For an example of a simple reference-counted type in Polymec, see the 
``st_func`` class in ``core/st_func.h`` (and its implementation in 
``core/st_func.c``). 

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
Experiment with these allocators to gain an understanding of their benefits, 
drawbacks, and general capabilities.

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

To formally document a type or a function, use Doxygen's markup:

http://www.doxygen.org/

In header files, describe your class types, structs, and enumerated types 
briefly and clearly. Build the Doxygen documentation to get an idea of what 
documentation typically looks like.

You don't need to put any documentation markup into implementation source 
files. Commenting your implementation code is always helpful, of course.

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
