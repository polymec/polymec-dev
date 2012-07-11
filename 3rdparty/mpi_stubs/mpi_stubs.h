# define MPI_COMM_WORLD 0

# define MPI_FAILURE 1
# define MPI_SUCCESS 0
# define MPI_ERR_IN_STATUS 2
# define MPI_ERR_TRUNCATE 3

# define MPI_SOURCE 1
# define MPI_TAG 2
# define MPI_COUNT 3

# define MPI_ANY_SOURCE -1
# define MPI_ANY_TAG -1

# define MPI_Comm int
# define MPI_Request int
# define MPI_Errhandler int

# define MPI_Datatype int
# define MPI_Op int

# define MPI_INT 1
# define MPI_FLOAT 2
# define MPI_DOUBLE 3
# define MPI_DOUBLE_PRECISION 3
# define MPI_BYTE 4
# define MPI_UNSIGNED 5
# define MPI_LONG 6
# define MPI_LONG_LONG 7
# define MPI_CHAR 8

# define MPI_SUM 1
# define MPI_MAX 2
# define MPI_MIN 3
# define MPI_PRODUCT 4

# define MPI_MAX_ERROR_STRING 1024
# define MPI_ERRORS_RETURN 0
# define MPI_ERRORS_ARE_FATAL 0

typedef struct 
{
  int MPI_ERROR;
} MPI_Status;

void MPI_Abort ( MPI_Comm comm, int ierror );
int MPI_Allgather ( void *sendbuf, int sendcount, MPI_Datatype sendtype,
  void *recvbuf, int recvcount, MPI_Datatype recvtype,
  MPI_Comm comm );
int MPI_Allgatherv ( void *sendbuf, int sendcount, MPI_Datatype sendtype,
  void *recvbuf, int *recvcounts, int *displs,
  MPI_Datatype recvtype, MPI_Comm comm );
int MPI_Allreduce ( void *sendbuf, void *recvbuf, int count,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm );
int MPI_Barrier ( MPI_Comm comm );
int MPI_Bcast ( void *data, int n, MPI_Datatype datatype, int node, 
  MPI_Comm comm );
int MPI_Cart_create ( MPI_Comm comm, int ndims, int dims[], int periods[],
  int reorder, MPI_Comm *comm_cart );
int MPI_Cart_get ( MPI_Comm comm, int ndims, int dims[], int periods[], 
  int coords[] );
int MPI_Cart_shift ( MPI_Comm comm, int dir, int disp, int *source, int *dest );
int MPI_Comm_dup ( MPI_Comm comm, MPI_Comm *comm_out );
int MPI_Comm_free ( MPI_Comm *comm );
int MPI_Comm_rank ( MPI_Comm comm, int *me );
int MPI_Comm_size ( MPI_Comm comm, int *nprocs );
int MPI_Comm_split ( MPI_Comm comm, int icolor, int ikey, MPI_Comm *new_comm );
int mpi_copy_byte ( char *data1, char *data2, int n );
int mpi_copy_double ( double *data1, double *data2, int n );
int mpi_copy_float ( float *data1, float *data2, int n );
int mpi_copy_int ( int *data1, int *data2, int n );
int MPI_Finalize ( void );
int MPI_Get_count ( MPI_Status status, MPI_Datatype datatype, int icount );
int MPI_Init ( int *argc, char **argv[] );
int MPI_Irecv ( void *buf, int count, MPI_Datatype datatype,
  int source, int tag, MPI_Comm comm, MPI_Request *request );
int MPI_Isend ( void *buf, int count, MPI_Datatype datatype,
  int dest, int tag, MPI_Comm comm, MPI_Request* request );
int MPI_Recv ( void *buf, int count, MPI_Datatype datatype,
  int source, int tag, MPI_Comm comm, MPI_Status *status );
int MPI_Reduce ( void *data1, void *data2, int n, MPI_Datatype datatype, 
  MPI_Op operation, int receiver, MPI_Comm comm );
int mpi_reduce_double ( double *data1, double *data2, int n, 
  MPI_Op operation );
int mpi_reduce_float ( float *data1, float *data2, int n, 
  MPI_Op operation );
int mpi_reduce_int ( int *data1, int *data2, int n, 
  MPI_Op operation );
int MPI_Reduce_scatter ( void *sendbuf, void *recvbuf, int recvcounts,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm );
int MPI_Rsend ( void *data, int n, MPI_Datatype datatype, int iproc, 
  int itag, MPI_Comm comm );
int MPI_Send ( void *buf, int count, MPI_Datatype datatype,
  int dest, int tag, MPI_Comm comm );
int MPI_Ssend ( void *buf, int count, MPI_Datatype datatype,
  int dest, int tag, MPI_Comm comm );
int MPI_Gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
               void *recvbuf, int recvcnt, MPI_Datatype recvtype, 
               int root, MPI_Comm comm);
int MPI_Gatherv(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
                void *recvbuf, int *recvcnts, int *displs, 
                MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Wait ( MPI_Request *request, MPI_Status *status );
int MPI_Waitall ( int icount, MPI_Request* requests, MPI_Status* statuses );
int MPI_Waitany ( int count, MPI_Request *request, int *index, MPI_Status *status );
int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status);
int MPI_Cancel(MPI_Request* request);
int MPI_Error_string(int errorcode, char *string, int *resultlen);
int MPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler *errhandler);
double MPI_Wtick ( void );
double MPI_Wtime ( void );
void timestamp ( void );

