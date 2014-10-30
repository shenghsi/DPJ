#include <mpi.h>
#include <stdio.h>
#include "MHD.h"
#include "global.h"
using namespace std;

int main(int argc, char* argv[]){
  double starttime, endtime;
  MPI_Init(&argc,&argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // number of order
  int no = 4;
  //  int ni,nj,nk;
  //  int nglobal[3] = {4,5,6};
  int nglobal[3] = {128,128,128};

  MHD mhd(rank, size, no,nglobal);
  
  mhd.initialize();
  
  while(mhd.t<10.){
  mhd.boundary();
  mhd.step();
  }
  mhd.UpdateGhostCell();

  // mhd.dump();
  // gather not working yet
  // mhd.gather();

  // mhd.hdfdump("mhd.hdf5");
  // mhd.hdfread("mhd.hdf5");
  // MPI_Barrier(MPI_COMM_WORLD);
  // starttime = MPI_Wtime();
  // //  mhd.hdfread("grid.h5");
  // mhd.hdfread("mhd1.hdf5");
  // cout<<rank<<" hdfread "<< MPI_Wtime() - starttime<<endl;

  // MPI_Barrier(MPI_COMM_WORLD);
  // starttime = MPI_Wtime();
  // mhd.hdfsave("mhd.hdf5");
  // cout<<rank<<" hdfsave "<< MPI_Wtime() - starttime<<endl;

  MPI_Barrier(MPI_COMM_WORLD);
  starttime = MPI_Wtime();
  mhd.hdfdump("mhd1.hdf5");
  cout<<rank<<" hdfdump "<< MPI_Wtime() - starttime<<endl;

  mhd.free_mpitype();
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();
  return 0;
}
