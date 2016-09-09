# Set up the compilers.

# If we're using MPI, install from source OpenMPI, using our nice 
# modern C compiler.
if [ $MPI -eq 1 ]; then
  wget https://www.open-mpi.org/software/ompi/v1.10/downloads/openmpi-1.10.4.tar.gz
  tar xzf openmpi-1.10.4.tar.gz
  pushd openmpi-1.10.4
  if [ "$CC" = "gcc" ]; then export CC="gcc-5"; fi
  if [ "$CC" = "clang" ]; then export CC="clang-3.6"; fi
  echo "Configuring OpenMPI..."
  env CC=$CC ./configure --disable-mpi-cxx --disable-java --prefix=/usr &> openmpi-config.log
  echo "Building and installing OpenMPI..."
  sudo make install &> openmpi-build.log
  popd
fi

# If we're building shared library and MPI support, install PETSc and HYPRE.
if [ $SHARED -eq 1 ] && [ $MPI -eq 1 ]; then
  . ./.travis/install-petsc-and-hypre.sh
fi 
