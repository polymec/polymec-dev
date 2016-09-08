# Set up the compilers.
if [ "$CC" = "gcc" ]; then export CC="gcc-5" CXX="g++-5"; fi
if [ "$CC" = "clang" ]; then export CC="clang-3.8" CXX="clang++-3.8"; fi

# If we're building shared library and MPI support, install PETSc and HYPRE.
if [ $SHARED -eq 1 ] && [ $MPI -eq 1 ]; then
  . ./.travis/install-petsc-and-hypre.sh
fi 
