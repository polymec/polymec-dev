# Install required software
brew update
brew install gcc@5 openmpi lcov
brew upgrade cmake

# Make sure the weird gfortran library links are in place.
ln -s /usr/local/lib/gcc/5/libgfortran.dylib /usr/local/lib/libgfortran.dylib
ln -s /usr/local/lib/gcc/5/libgfortran.a /usr/local/lib/libgfortran.a

# If we're building shared library support, install PETSc and HYPRE.
if [ $SHARED -eq 1 ]; then
  . ./.travis/install-petsc-and-hypre.sh
fi 
