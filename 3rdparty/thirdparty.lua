-- This metadata describes the various third-party libraries in polymec,
-- and is exposed in the thirdparty module.

thirdparty = {

cmocka = {
  version = '1.1.3',
  copyright = 'Copyright (c) 2008 Google Inc.\n' ..
              'Copyright (c) 2013-2015,2016,2018, Andreas Schneider',
  license = 'Apache',
  url = 'https://cmocka.org',
  date_retrieved = {month=11, day=12, year=2018}
},

hdf5 = {
  version = '1.10.4',
  copyright = 'Copyright (c) 2006-2018, The HDF Group.\n' .. 
              'Copyright (c) 1998-2006, The Board of Trustees of the University of Illinois.\n' ..
              'All rights reserved.',
  url = 'https://www.hdfgroup.org/downloads/hdf5/source-code',
  license = 'BSD',
  date_retrieved = {month=11, day=12, year=2018}
},

jrs_predicates = {
  version = '05-18-1996',
  url = 'http://www.cs.cmu.edu/~quake/robust.html',
  license = 'public domain',
  date_retrieved = {month=5, day=18, year=1996},
  notes = 'This is predicates.c, released to the public domain by Jonathan Richard ' ..
          'Shewchuk (as explained in the source). The only modification made to this ' ..
          'file for inclusion in polymec was to rename compress() to compress_expansion() ' .. 
          'to avoid symbol interference with zlib.'
},

libarena = {
  version = '0.3.7',
  copyright = 'Copyright (c) 2006  William Ahern',
  url = 'http://www.25thandclement.com/~william/projects/libarena.html',
  license = 'BSD',
  date_retrieved = {month=1, day=20, year=2014}
},

linenoise = {
  version = '1.0',
  copyright = 'Copyright (c) 2010-2014, Salvatore Sanfilippo\n' .. 
              'Copyright (c) 2010-2013, Pieter Noordhuis\n' .. 
              'All rights reserved.',
  license = 'BSD',
  url = 'https://github.com/antirez/linenoise',
  date_retrieved = {month=4, day=13, year=2015},
  notes = '# Modified (see commit 56baf39)'
},

lua = {
  version = '5.3.5',
  copyright = 'Copyright (c) 1994-2017 Lua.org, PUC-Rio.',
  url = 'http://www.lua.org',
  license = 'MIT',
  date_retrieved = {month=7, day=25, year=2018}
},

scotch = {
  version = '6.0.6',
  copyright = 'Copyright (c) 2004,2007,2008,2010-2012,2014,2018 IPB, Universite de Bordeaux, IN',
  url = 'http://gforge.inria.fr/projects/scotch',
  license = 'CeCILL-C',
  date_retrieved = {month=7, day=25, year=2018}
},

usilo = {
  version = '5.0',
  copyright = 'Copyright (c) 1994-2010, Lawrence Livermore National Security, LLC.\n' .. 
              'LLNL-CODE-425250.\n' .. 
              'All rights reserved.',
  license = 'BSD',
  date_retrieved = {month=4, day=20, year=2018},
  notes = 'Forked from Silo 4.10.4, a limited release from ' .. 
          'https://wci.llnl.gov/simulation/computer-codes/silo. ' .. 
          'Silo is abandonware, so "microsilo" is a fork intended ' .. 
          'only for use with polymec.'
},

sundials = {
  version = '4.0.1',
  copyright = 'Copyright (c) 2002-2016, Lawrence Livermore National Security.\n' .. 
              'Produced at the Lawrence Livermore National Laboratory.\n' .. 
              'Written by A.C. Hindmarsh, D.R. Reynolds, R. Serban, C.S. Woodward,\n' .. 
              'S.D. Cohen, A.G. Taylor, S. Peles, L.E. Banks, and D. Shumaker.\n' ..
              'LLNL-CODE-667205    (ARKODE)\n' .. 
              'UCRL-CODE-155951    (CVODE)\n' .. 
              'UCRL-CODE-155950    (CVODES)\n' .. 
              'UCRL-CODE-155952    (IDA)\n' .. 
              'UCRL-CODE-237203    (IDAS)\n' .. 
              'LLNL-CODE-665877    (KINSOL)\n' .. 
              'All rights reserved.',
  license = 'BSD',
  url = 'https://computation.llnl.gov/projects/sundials/sundials-software',
  date_retrieved = {month=12, day=23, year=2018},
  citation = 'Alan C. Hindmarsh, Peter N. Brown, Keith E. Grant, Steven L. Lee, Radu ' .. 
             'Serban, Dan E. Shumaker, and Carol S. Woodward. 2005. SUNDIALS: Suite of ' .. 
             'nonlinear and differential/algebraic equation solvers. ACM Trans. Math. Softw. ' .. 
             '31, 3 (September 2005), 363-396. DOI=http://dx.doi.org/10.1145/1089014.1089020'
},

zlib = {
  version = '1.2.11',
  copyright = 'Copyright (c) 1995-2017 Jean-loup Gailly and Mark Adler',
  license = 'zlib',
  url = 'http://www.zlib.net',
  date_retrieved = {month=9, day=21, year=2018}
}

}
