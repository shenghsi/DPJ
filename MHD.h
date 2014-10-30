#include <mpi.h>
class MHD{
 public:
  MHD(int, int, int, int *);
  ~MHD();
  void initialize();
  void boundary();
  void step();
  void dump();
  void free_mpitype();
  void UpdateGhostCell();
  void gather();
  void hdfdump(const char*);
  void hdfread(const char*);
  void hdfsave(const char*);
  void hdfload(const char*);
  float t;
  
 private:
  //MPI
  int rank, size;
  MPI_Comm comm_cart;
  int nbr_left[3],nbr_right[3];
  int coord[3];

  int proc_ndims;
  int proc_dims[3];

  
  // Datatype for variables at cell centers of ghost cells
  MPI_Datatype xall_send_left;
  MPI_Datatype xall_recv_right;
  MPI_Datatype xall_send_right;
  MPI_Datatype xall_recv_left;
  MPI_Datatype yall_send_left;
  MPI_Datatype yall_recv_right;
  MPI_Datatype yall_send_right;
  MPI_Datatype yall_recv_left;
  MPI_Datatype zall_send_left;
  MPI_Datatype zall_recv_right;
  MPI_Datatype zall_send_right;
  MPI_Datatype zall_recv_left;

  //order of scheme
  int NO;
  // grid size
  int ni_global,nj_global,nk_global;
  int NO2; // number of ghost cells
  //local sizes:
  int NI, NJ, NK;  //local active cell size
  int NIT, NJT, NKT; // total cell centers
  int NITf, NJTf, NKTf; // total cell faces
  //bound
  int ilobound;
  int jlobound;
  int klobound;
  int ihibound;
  int jhibound;
  int khibound;

  // index
  int ial; // first active cell
  int iah; // last active cell
  int jal; // first active cell
  int jah; // last active cell
  int kal;
  int kah;

  // offset for the first active cell of each processor respect to the first active cell of global array
  int offset_global[3];

  // grid
  float ***x, ***y, ***z;
  // cell center variables
  float ***rho,***p;
  float ***vx,***vy,***vz;
  float ***bx,***by,***bz;
  //b on face
  float ***bi, ***bj, ***bk;
  //e on edge
  //  float ***ei, ***ej, ***ek;
  // conserved variables:
  //  float ***rhovx,***rhovy,***rhovz,***eng;

  void courant();
  void reconstruct();
  void hydro_flux();
  void magnetic_flux();
  void get_conserved_variables();
};
