# Install necessary software.
sudo apt-get update -qq
sudo apt-get install -y cmake gcc clang libopenmpi-dev openmpi-bin liblapack-dev gfortran git valgrind

# If we're building shared library and MPI support, install PETSc and HYPRE.
if [ $SHARED -eq 1 ] && [ $MPI -eq 1 ]; then
  . ./.travis/install-petsc-and-hypre.sh
fi 
