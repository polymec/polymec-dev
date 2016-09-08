# Set up the compilers.
# If we're building shared library and MPI support, install PETSc and HYPRE.
if [ $SHARED -eq 1 ] && [ $MPI -eq 1 ]; then
  . ./.travis/install-petsc-and-hypre.sh
fi 
