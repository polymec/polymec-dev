# Go get PETSc 3.6.x and build it.
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.6.3.tar.gz
tar xzf petsc-lite-3.6.3.tar.gz
pushd $PETSC_DIR
./configure --with-mpi=$MPI --with-debugging=$DEBUG --with-shared-libraries=1 --with-64-bit-indices=1
make
ln -s $PETSC_DIR/lib/petsc/conf $PETSC_DIR/conf
ln -s $PETSC_DIR/include/petsc/finclude $PETSC_DIR/include/finclude
popd

# Do the same for HYPRE.
wget http://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/hypre-2.10.1.tar.gz
tar xzf hypre-2.10.1.tar.gz
pushd $HYPRE_DIR
if [ $MPI -eq 0 ]; then
  HYPRE_SEQ=ON
else
  HYPRE_SEQ=OFF
fi 
cmake .. -DHYPRE_SHARED=ON -DHYPRE_BIGINT=ON -DHYPRE_SEQUENTIAL=$HYPRE_SEQ -DHYPRE_USING_FEI=OFF
make -j4
popd

