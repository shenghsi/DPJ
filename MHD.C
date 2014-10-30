#include "MHD.h"
#include <math.h>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include "arrays.h"
using namespace std;

//#include <stdio.h>
MHD::MHD(int nrank, int nsize, int no, int *global_size)
 :NO(no)
{
  ni_global = global_size[0];
  nj_global = global_size[1];
  nk_global = global_size[2];
  //
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // MPI_Comm_size(MPI_COMM_WORLD, &size);
  rank = nrank;
  size = nsize;
  // Cartesian processors layout
  // number of dimensions for domain decomposition:
  proc_ndims = 3;
  // number of procs in each direction:
  proc_dims[0] = 4;
  proc_dims[1] = 4;
  proc_dims[2] = 2;

  if(size != proc_dims[0]*proc_dims[1]*proc_dims[2]){
    cout<<"proc number does not agree with internal setup"<<endl;
    exit(1);
  }

  // Cartesian layout for procs:
  int periods[3] = {0,0,0};
  int reorder = 0;
  // New communicator comm_cart with processes ordered in a Cartesian grid
  MPI_Cart_create(MPI_COMM_WORLD, proc_ndims, proc_dims, periods, reorder, &comm_cart);
  // processor coordinates coord
  MPI_Cart_coords(comm_cart, rank, 3, coord);
  // find neighbors for each processor; MPI_PROC_NULL is set if no nerghbors are found
  int dim;
  dim = 0;
  MPI_Cart_shift(comm_cart,dim,1,&nbr_left[dim],&nbr_right[dim]);
  dim = 1;
  MPI_Cart_shift(comm_cart,dim,1,&nbr_left[dim],&nbr_right[dim]);
  dim = 2;
  MPI_Cart_shift(comm_cart,dim,1,&nbr_left[dim],&nbr_right[dim]);
  //  cout<<"rank= "<<rank<<" dim0: "<<nbr_left[0]<<" "<<nbr_right[0]<<" dim1: "<<nbr_left[1]<<" "<<nbr_right[1]<<" dim2:"<<nbr_left[2]<<" "<<nbr_right[2]<<endl;

  // figure out if the processor is on any boundary
  ilobound = 0;
  jlobound = 0;
  klobound = 0;
  ihibound = 0;
  jhibound = 0;
  khibound = 0;
  if(nbr_left[0] == MPI_PROC_NULL)
    ilobound = 1;
  if(nbr_left[1] == MPI_PROC_NULL)
    jlobound = 1;
  if(nbr_left[2] == MPI_PROC_NULL)
    klobound = 1;
  if(nbr_right[0] == MPI_PROC_NULL)
    ihibound = 1;
  if(nbr_right[1] == MPI_PROC_NULL)
    jhibound = 1;
  if(nbr_right[2] == MPI_PROC_NULL)
    khibound = 1;

  // cout<<"rank = " << rank<<" "
  //     <<"coord = "<<coord[0]<<" "<<coord[1]<<" "<<coord[2]<<" "
  //     <<"lobound = "<<ilobound<<" "<<jlobound<<" "<<klobound<<" "
  //     <<"hibound = "<<ihibound<<" "<<jhibound<<" "<<khibound<<" "
  //     <<endl;

  // find local size and offsets for the first active cell
  int local_size[3];
  for(int i=0;i<3;i++){
    local_size[i] = global_size[i]/proc_dims[i];
    int extra;
    extra = (global_size[i] - proc_dims[i]*local_size[i]);
    if(coord[i] < extra) {
      local_size[i]++;
      offset_global[i] = coord[i]*local_size[i];
    }
    else {
      offset_global[i] = coord[i]*local_size[i] + extra;
    }
  }
  // cout<<" rank = "<<rank
  //     <<" coord = "<<coord[0]<<" "<<coord[1]<<" "<<coord[2]
  //     <<" local_size = "<< local_size[0]<<" "<< local_size[1]<<" "<< local_size[2]
  //     <<" offset_global = " << offset_global[0]<<" "<< offset_global[1]<<" "<< offset_global[2]
  //     <<endl;

  //
  t=0.;
  NI = local_size[0];
  NJ = local_size[1];
  NK = local_size[2];

  NO2 = NO/2;
  NIT = NI + NO;
  NJT = NJ + NO;
  NKT = NK + NO;

  NITf = NIT+1;
  NJTf = NJT+1;
  NKTf = NKT+1;

  ial = NO2;
  iah = NI+NO2-1;
  jal = NO2;
  jah = NJ+NO2-1;
  kal = NO2;
  kah = NK+NO2-1;

  x = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));
  y = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));
  z = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));
  rho = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));
  p = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));
  vx = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));
  vy = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));
  vz = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));
  bx = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));
  by = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));
  bz = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));

  bi = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));
  bj = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));
  bk = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));
  // ei = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));
  // ej = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));
  // ek = (float***) New3DArray(NITf,NJTf,NKTf,sizeof(float));


  // //**************************
  for(int i=0; i<NITf; i++)
    for(int j=0; j<NJTf; j++)
      for(int k=0; k<NKTf; k++){
  	rho[i][j][k] = rank*1000+ i*100+j*10+k;
  	//	p[i][j][k] = rank*1000+ i*100+j*10+k;
  	p[i][j][k] =  (offset_global[0]+i-NO2)*10000 + (offset_global[1]+j-NO2)*100 + offset_global[2]+k-NO2;
      }

  int tag = 1;
  MPI_Status status;

  int local_size_with_ghost[3] = {NITf,NJTf,NKTf};

  //x direction:
  //x centers
  int x_sub_size[3] = {NO2,NJ,NK};
  //send to Left / recv from Right
  int xstart_send_left[3]  = {ial,  jal,kal};
  int xstart_recv_right[3] = {iah+1,jal,kal};
  //send to Right / recv from Left
  int xstart_send_right[3] = {iah+1-NO2,jal,kal};
  int xstart_recv_left[3]  = {0,        jal,kal};
  MPI_Datatype xtype_send_left;
  MPI_Datatype xtype_recv_right;
  MPI_Datatype xtype_send_right;
  MPI_Datatype xtype_recv_left;
  MPI_Type_create_subarray(3, local_size_with_ghost, x_sub_size, xstart_send_left, MPI_ORDER_C,MPI_FLOAT, &xtype_send_left);
  MPI_Type_commit(&xtype_send_left);
  MPI_Type_create_subarray(3, local_size_with_ghost, x_sub_size, xstart_recv_right, MPI_ORDER_C, MPI_FLOAT, &xtype_recv_right);
  MPI_Type_commit(&xtype_recv_right);
  MPI_Type_create_subarray(3, local_size_with_ghost, x_sub_size, xstart_send_right, MPI_ORDER_C,MPI_FLOAT, &xtype_send_right);
  MPI_Type_commit(&xtype_send_right);
  MPI_Type_create_subarray(3, local_size_with_ghost, x_sub_size, xstart_recv_left, MPI_ORDER_C, MPI_FLOAT, &xtype_recv_left);
  MPI_Type_commit(&xtype_recv_left);
  //x face variables: bi, ej, ek
  //send to Left / recv from Right
  xstart_send_left[0]  = ial+1;
  xstart_recv_right[0] = iah+2;
  //send to Right / recv from Left
  //  xstart_send_right[0] = iah+1-NO2;
  //  xstart_recv_left[3]  = {0,        jal,kal};
  MPI_Datatype xface_send_left;
  MPI_Datatype xface_recv_right;
  MPI_Datatype xface_send_right;
  MPI_Datatype xface_recv_left;
  MPI_Type_create_subarray(3, local_size_with_ghost, x_sub_size, xstart_send_left, MPI_ORDER_C,MPI_FLOAT, &xface_send_left);
  MPI_Type_commit(&xface_send_left);
  MPI_Type_create_subarray(3, local_size_with_ghost, x_sub_size, xstart_recv_right, MPI_ORDER_C, MPI_FLOAT, &xface_recv_right);
  MPI_Type_commit(&xface_recv_right);
  MPI_Type_create_subarray(3, local_size_with_ghost, x_sub_size, xstart_send_right, MPI_ORDER_C,MPI_FLOAT, &xface_send_right);
  MPI_Type_commit(&xface_send_right);
  MPI_Type_create_subarray(3, local_size_with_ghost, x_sub_size, xstart_recv_left, MPI_ORDER_C, MPI_FLOAT, &xface_recv_left);
  MPI_Type_commit(&xface_recv_left); 

  //y direction:
  int y_sub_size[3] = {NI,NO2,NK};
  //send to Left / recv from Right
  int ystart_send_left[3]  = {ial,jal,  kal};
  int ystart_recv_right[3] = {ial,jah+1,kal};
  //send to Right / recv from Left
  int ystart_send_right[3] = {ial,jah+1-NO2,kal};
  int ystart_recv_left[3]  = {ial,0,        kal};
  MPI_Datatype ytype_send_left;
  MPI_Datatype ytype_recv_right;
  MPI_Datatype ytype_send_right;
  MPI_Datatype ytype_recv_left;
  MPI_Type_create_subarray(3, local_size_with_ghost, y_sub_size, ystart_send_left, MPI_ORDER_C,MPI_FLOAT, &ytype_send_left);
  MPI_Type_commit(&ytype_send_left);
  MPI_Type_create_subarray(3, local_size_with_ghost, y_sub_size, ystart_recv_right, MPI_ORDER_C, MPI_FLOAT, &ytype_recv_right);
  MPI_Type_commit(&ytype_recv_right);
  MPI_Type_create_subarray(3, local_size_with_ghost, y_sub_size, ystart_send_right, MPI_ORDER_C,MPI_FLOAT, &ytype_send_right);
  MPI_Type_commit(&ytype_send_right);
  MPI_Type_create_subarray(3, local_size_with_ghost, y_sub_size, ystart_recv_left, MPI_ORDER_C, MPI_FLOAT, &ytype_recv_left);
  MPI_Type_commit(&ytype_recv_left);
  //y face: bj, ek, ei
  //send to Left / recv from Right
  ystart_send_left[1]  = jal+1;
  ystart_recv_right[1] = jah+2;
  //send to Right / recv from Left
  //  ystart_send_right[3] = {ial,jah+1-NO2,kal};
  //  ystart_recv_left[3]  = {ial,0,        kal};
  MPI_Datatype yface_send_left;
  MPI_Datatype yface_recv_right;
  MPI_Datatype yface_send_right;
  MPI_Datatype yface_recv_left;
  MPI_Type_create_subarray(3, local_size_with_ghost, y_sub_size, ystart_send_left, MPI_ORDER_C,MPI_FLOAT, &yface_send_left);
  MPI_Type_commit(&yface_send_left);
  MPI_Type_create_subarray(3, local_size_with_ghost, y_sub_size, ystart_recv_right, MPI_ORDER_C, MPI_FLOAT, &yface_recv_right);
  MPI_Type_commit(&yface_recv_right);
  MPI_Type_create_subarray(3, local_size_with_ghost, y_sub_size, ystart_send_right, MPI_ORDER_C,MPI_FLOAT, &yface_send_right);
  MPI_Type_commit(&yface_send_right);
  MPI_Type_create_subarray(3, local_size_with_ghost, y_sub_size, ystart_recv_left, MPI_ORDER_C, MPI_FLOAT, &yface_recv_left);
  MPI_Type_commit(&yface_recv_left);

  //z direction:
  int z_sub_size[3] = {NI,NJ,NO2};
  //send to Left / recv from Right
  int zstart_send_left[3]  = {ial,jal,kal};
  int zstart_recv_right[3] = {ial,jal,kah+1};
  //send to Right / recv from Left
  int zstart_send_right[3] = {ial,jal,kah+1-NO2};
  int zstart_recv_left[3]  = {ial,jal,0};
  MPI_Datatype ztype_send_left;
  MPI_Datatype ztype_recv_right;
  MPI_Datatype ztype_send_right;
  MPI_Datatype ztype_recv_left;
  MPI_Type_create_subarray(3, local_size_with_ghost, z_sub_size, zstart_send_left, MPI_ORDER_C,MPI_FLOAT, &ztype_send_left);
  MPI_Type_commit(&ztype_send_left);
  MPI_Type_create_subarray(3, local_size_with_ghost, z_sub_size, zstart_recv_right, MPI_ORDER_C, MPI_FLOAT, &ztype_recv_right);
  MPI_Type_commit(&ztype_recv_right);
  MPI_Type_create_subarray(3, local_size_with_ghost, z_sub_size, zstart_send_right, MPI_ORDER_C,MPI_FLOAT, &ztype_send_right);
  MPI_Type_commit(&ztype_send_right);
  MPI_Type_create_subarray(3, local_size_with_ghost, z_sub_size, zstart_recv_left, MPI_ORDER_C, MPI_FLOAT, &ztype_recv_left);
  MPI_Type_commit(&ztype_recv_left);
  //z face: bk, ei, ej
  //send to Left / recv from Right
  zstart_send_left[2]  = kal+1;
  zstart_recv_right[2] = kah+2;
  //send to Right / recv from Left
  //  zstart_send_right[3] = {ial,jal,kah+1-NO2};
  //  zstart_recv_left[3]  = {ial,jal,0};
  MPI_Datatype zface_send_left;
  MPI_Datatype zface_recv_right;
  MPI_Datatype zface_send_right;
  MPI_Datatype zface_recv_left;
  MPI_Type_create_subarray(3, local_size_with_ghost, z_sub_size, zstart_send_left, MPI_ORDER_C,MPI_FLOAT, &zface_send_left);
  MPI_Type_commit(&zface_send_left);
  MPI_Type_create_subarray(3, local_size_with_ghost, z_sub_size, zstart_recv_right, MPI_ORDER_C, MPI_FLOAT, &zface_recv_right);
  MPI_Type_commit(&zface_recv_right);
  MPI_Type_create_subarray(3, local_size_with_ghost, z_sub_size, zstart_send_right, MPI_ORDER_C,MPI_FLOAT, &zface_send_right);
  MPI_Type_commit(&zface_send_right);
  MPI_Type_create_subarray(3, local_size_with_ghost, z_sub_size, zstart_recv_left, MPI_ORDER_C, MPI_FLOAT, &zface_recv_left);
  MPI_Type_commit(&zface_recv_left);

  // MPI_Sendrecv(&rho[0][0][0], 1, xtype_send_left, nbr_left[0], tag, &rho[0][0][0], 1, xtype_recv_right, nbr_right[0], MPI_ANY_TAG, comm_cart, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(&rho[0][0][0], 1, xtype_send_right, nbr_right[0], tag, &rho[0][0][0], 1, xtype_recv_left, nbr_left[0], MPI_ANY_TAG, comm_cart, MPI_STATUS_IGNORE);
  // if(rank==1)
  //   MPI_Send(&rho[0][0][0], 1, xtype_send_left, 0, tag, comm_cart);
  // else if(rank==0)
  //   MPI_Recv(&rho[0][0][0], 1, xtype_recv_right, 1, MPI_ANY_TAG, comm_cart, &status);
  // MPI_Send(&rho[0][0][0], 1, xtype_send_left, nbr_left[0], tag, comm_cart);
  // MPI_Recv(&rho[0][0][0], 1, xtype_recv_right, nbr_right[0], MPI_ANY_TAG, comm_cart, &status);

  // Create data type including all the variables that need exchange:
  int nvar = 10;
  int blen[nvar];
  MPI_Aint disp[nvar];

  for(int n=0;n<nvar;n++)
    blen[n] = 1;

  // get address for different variables and calculate the relative address displacements:
  //***************************************************
  // x direction:
  MPI_Get_address(&rho[0][0][0],&disp[0]);
  MPI_Get_address(&p[0][0][0],&disp[1]);
  MPI_Get_address(&vx[0][0][0],&disp[2]);
  MPI_Get_address(&vy[0][0][0],&disp[3]);
  MPI_Get_address(&vz[0][0][0],&disp[4]);
  MPI_Get_address(&bx[0][0][0],&disp[5]);
  MPI_Get_address(&by[0][0][0],&disp[6]);
  MPI_Get_address(&bz[0][0][0],&disp[7]);
  //only bj, bk needed for x direction
  MPI_Get_address(&bj[0][0][0],&disp[8]);
  MPI_Get_address(&bk[0][0][0],&disp[9]);
  for(int n=1;n<nvar;n++){
    disp[n] = disp[n] - disp[0];
  }
  disp[0] = 0;

  //MPI Datatype for variables at cell centers in x direction:
  MPI_Type_hindexed(nvar, blen, disp, xtype_send_left, &xall_send_left);
  MPI_Type_hindexed(nvar, blen, disp, xtype_recv_right, &xall_recv_right);
  MPI_Type_commit(&xall_send_left);
  MPI_Type_commit(&xall_recv_right);

  MPI_Type_hindexed(nvar, blen, disp, xtype_send_right, &xall_send_right);
  MPI_Type_hindexed(nvar, blen, disp, xtype_recv_left, &xall_recv_left);
  MPI_Type_commit(&xall_send_right);
  MPI_Type_commit(&xall_recv_left);

  //***************************************************
  // y direction:
  MPI_Get_address(&rho[0][0][0],&disp[0]);
  MPI_Get_address(&p[0][0][0],&disp[1]);
  MPI_Get_address(&vx[0][0][0],&disp[2]);
  MPI_Get_address(&vy[0][0][0],&disp[3]);
  MPI_Get_address(&vz[0][0][0],&disp[4]);
  MPI_Get_address(&bx[0][0][0],&disp[5]);
  MPI_Get_address(&by[0][0][0],&disp[6]);
  MPI_Get_address(&bz[0][0][0],&disp[7]);
  //only bk, bi needed for y direction
  MPI_Get_address(&bk[0][0][0],&disp[8]);
  MPI_Get_address(&bi[0][0][0],&disp[9]);
  for(int n=1;n<nvar;n++){
    disp[n] = disp[n] - disp[0];
  }
  disp[0] = 0;
  //MPI Datatype for variables at cell centers in y direction:
  MPI_Type_hindexed(nvar, blen, disp, ytype_send_left, &yall_send_left);
  MPI_Type_hindexed(nvar, blen, disp, ytype_recv_right, &yall_recv_right);
  MPI_Type_commit(&yall_send_left);
  MPI_Type_commit(&yall_recv_right);

  MPI_Type_hindexed(nvar, blen, disp, ytype_send_right, &yall_send_right);
  MPI_Type_hindexed(nvar, blen, disp, ytype_recv_left, &yall_recv_left);
  MPI_Type_commit(&yall_send_right);
  MPI_Type_commit(&yall_recv_left);

  //***************************************************
  // z direction:
  MPI_Get_address(&rho[0][0][0],&disp[0]);
  MPI_Get_address(&p[0][0][0],&disp[1]);
  MPI_Get_address(&vx[0][0][0],&disp[2]);
  MPI_Get_address(&vy[0][0][0],&disp[3]);
  MPI_Get_address(&vz[0][0][0],&disp[4]);
  MPI_Get_address(&bx[0][0][0],&disp[5]);
  MPI_Get_address(&by[0][0][0],&disp[6]);
  MPI_Get_address(&bz[0][0][0],&disp[7]);
  //only bi, bj needed for z direction
  MPI_Get_address(&bi[0][0][0],&disp[8]);
  MPI_Get_address(&bj[0][0][0],&disp[9]);
  for(int n=1;n<nvar;n++){
    disp[n] = disp[n] - disp[0];
  }
  disp[0] = 0;
  //MPI Datatype for variables at cell centers in z direction:  
  MPI_Type_hindexed(nvar, blen, disp, ztype_send_left, &zall_send_left);
  MPI_Type_hindexed(nvar, blen, disp, ztype_recv_right, &zall_recv_right);
  MPI_Type_commit(&zall_send_left);
  MPI_Type_commit(&zall_recv_right);

  MPI_Type_hindexed(nvar, blen, disp, ztype_send_right, &zall_send_right);
  MPI_Type_hindexed(nvar, blen, disp, ztype_recv_left, &zall_recv_left);
  MPI_Type_commit(&zall_send_right);
  MPI_Type_commit(&zall_recv_left);

  // test exchange
  // MPI_Sendrecv(&rho[0][0][0], 1, xall_send_left, nbr_left[0], tag, &rho[0][0][0], 1, xall_recv_right, nbr_right[0], MPI_ANY_TAG, comm_cart, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(&rho[0][0][0], 1, xall_send_right, nbr_right[0], tag, &rho[0][0][0], 1, xall_recv_left, nbr_left[0], MPI_ANY_TAG, comm_cart, MPI_STATUS_IGNORE);

  // MPI_Sendrecv(&rho[0][0][0], 1, yall_send_left, nbr_left[1], tag, &rho[0][0][0], 1, yall_recv_right, nbr_right[1], MPI_ANY_TAG, comm_cart, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(&rho[0][0][0], 1, yall_send_right, nbr_right[1], tag, &rho[0][0][0], 1, yall_recv_left, nbr_left[1], MPI_ANY_TAG, comm_cart, MPI_STATUS_IGNORE);

  // MPI_Sendrecv(&rho[0][0][0], 1, zall_send_left, nbr_left[2], tag, &rho[0][0][0], 1, zall_recv_right, nbr_right[2], MPI_ANY_TAG, comm_cart, MPI_STATUS_IGNORE);
  // MPI_Sendrecv(&rho[0][0][0], 1, zall_send_right, nbr_right[2], tag, &rho[0][0][0], 1, zall_recv_left, nbr_left[2], MPI_ANY_TAG, comm_cart, MPI_STATUS_IGNORE);

  //x:
  MPI_Type_free(&xtype_send_left);  
  MPI_Type_free(&xtype_recv_right);  
  MPI_Type_free(&xtype_send_right);  
  MPI_Type_free(&xtype_recv_left);  
  //y:
  MPI_Type_free(&ytype_send_left);  
  MPI_Type_free(&ytype_recv_right);  
  MPI_Type_free(&ytype_send_right);  
  MPI_Type_free(&ytype_recv_left);  
  //z:
  MPI_Type_free(&ztype_send_left);  
  MPI_Type_free(&ztype_recv_right);  
  MPI_Type_free(&ztype_send_right);  
  MPI_Type_free(&ztype_recv_left);  

  //face
  //x:
  MPI_Type_free(&xface_send_left);  
  MPI_Type_free(&xface_recv_right);  
  MPI_Type_free(&xface_send_right);  
  MPI_Type_free(&xface_recv_left);  
  //y:
  MPI_Type_free(&yface_send_left);  
  MPI_Type_free(&yface_recv_right);  
  MPI_Type_free(&yface_send_right);  
  MPI_Type_free(&yface_recv_left);  
  //z:
  MPI_Type_free(&zface_send_left);  
  MPI_Type_free(&zface_recv_right);  
  MPI_Type_free(&zface_send_right);  
  MPI_Type_free(&zface_recv_left);  


  // if(rank==0)
  //   for(int k=0; k<NKTf; k++){
  //     for(int i=0; i<NITf; i++){
  // 	for(int j=0; j<NJTf; j++){
  // 	  cout<<p[i][j][k]<<" ";
  // 	}
  // 	cout<<endl;
  //     }
  //     cout<<endl;
  //   }

  //====================================
  // use MPI_Type_indexed
  //====================================
  // int buf_send,buf_recv;
  // buf_recv = -999;
  // buf_send = rank;

  // MPI_Send(&buf_send,1,MPI_INT, nbr_left[0], sendtag, comm_cart);
  // MPI_Recv(&buf_recv,1,MPI_INT, nbr_right[0], recvtag, comm_cart, status);

  // cout<<"rank "<< rank <<" buf_send " <<buf_send<<" buf_recv "<<buf_recv<<endl;
  // MPI_Type_free(&x_recv_backward);  
  
  // MPI_Datatype type_x;
  // int blocklen[NO2*NJ], disp[NO2*NJ];
  // for(int i=0;i<NO2;i++){
  //   for(int j=0;j<NJ;j++){
  //     int n = i*NJ + j;
  //     blocklen[n] = NK;
  //     disp[n] = i*NJTf*NKTf + j*NKTf;
  //   }
  // }
  // MPI_Type_indexed(NO2*NJ, blocklen, disp, MPI_FLOAT, &type_x);
  // MPI_Type_commit(&type_x);
  // int tag;
  // MPI_Send(&rho[ial][jal][kal], 1, type_x, nbr_left[0], tag, comm_cart);
  // MPI_Status status;
  // MPI_Recv(&rho[0][jal][kal], 1, type_x, nbr_right[0], tag, comm_cart, &status);
  // MPI_Type_free(&type_x);

}

MHD::~MHD(){
  Delete3DArray((void***)x);
  Delete3DArray((void***)rho);
  Delete3DArray((void***)p);
  Delete3DArray((void***)vx);
  Delete3DArray((void***)vy);
  Delete3DArray((void***)vz);
  Delete3DArray((void***)bx);
  Delete3DArray((void***)by);
  Delete3DArray((void***)bz);
  Delete3DArray((void***)bi);
  Delete3DArray((void***)bj);
  Delete3DArray((void***)bk);
  // Delete3DArray((void***)ei);
  // Delete3DArray((void***)ej);
  // Delete3DArray((void***)ek);
}
void MHD::free_mpitype(){
  MPI_Type_free(&xall_send_left);
  MPI_Type_free(&xall_recv_right);
  MPI_Type_free(&xall_send_right);
  MPI_Type_free(&xall_recv_left);
  MPI_Type_free(&yall_send_left);
  MPI_Type_free(&yall_recv_right);
  MPI_Type_free(&yall_send_right);
  MPI_Type_free(&yall_recv_left);
  MPI_Type_free(&zall_send_left);
  MPI_Type_free(&zall_recv_right);
  MPI_Type_free(&zall_send_right);
  MPI_Type_free(&zall_recv_left);
}
void MHD::initialize(){

  for(int i=0; i<NITf; i++)
    for(int j=0; j<NJTf; j++)
      for(int k=0; k<NKTf; k++){
  	x[i][j][k] = (offset_global[0]+i-NO2)*(offset_global[0]+i-NO2) + float(offset_global[1]+j-NO2) + sqrt(float(offset_global[2]+k-NO2));
  	y[i][j][k] = rank*1000+ i*100+j*10+k;
  	z[i][j][k] =  (offset_global[0]+i-NO2)*10000 + (offset_global[1]+j-NO2)*100 + offset_global[2]+k-NO2;
      }


  // float dx;
  // dx = (Lx[1]-Lx[0])/NI;
  // x[0] = Lx[0] - NO2*dx;
  // for(int i=1; i<NITf; i++)
  //   x[i] = x[i-1] + dx;
}

void MHD::boundary(){

}

void MHD::step(){
  t=t+0.1;
  courant();
  reconstruct();
  hydro_flux();
  magnetic_flux();
  get_conserved_variables();
}

void MHD::dump(){
  clock_t begin, end;
  double time_spent;


  if(rank==0)
    for(int i=0; i<NITf; i++){
      for(int j=0; j<NJTf; j++){
	for(int k=0; k<NKTf; k++){
  	  cout<<p[i][j][k]<<" ";
	}
  	cout<<endl;
      }
      cout<<endl;
    }


  // begin = clock();
  // /* here, do your time-consuming job */
  // for(int i=0;i<10;i++)
  //   cout<< erf((float) i) << endl;
  //   //    erf(t);
  // end = clock();
  // time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  // //  double x = 0.9;
  // cout<<"time_spent"<<setprecision(10)<< time_spent << endl;

  // for(int i=0;i<NITf;i++)
  //   cout<<x[i][0][0]<<" ";
  // cout<<endl;

}

void MHD::courant(){
}
void MHD::reconstruct(){
}
void MHD::hydro_flux(){
}
void MHD::magnetic_flux(){
}
void MHD::get_conserved_variables(){
}

// rho, p, vx,vy,vz, bx,by,bz ghost cells are filled for all 3 directions
// bi,bj,bk ghost cells are filled only for 2 directions
void MHD::UpdateGhostCell(){
  int tag;

  MPI_Sendrecv(&rho[0][0][0], 1, xall_send_left, nbr_left[0], tag, &rho[0][0][0], 1, xall_recv_right, nbr_right[0], MPI_ANY_TAG, comm_cart, MPI_STATUS_IGNORE);
  MPI_Sendrecv(&rho[0][0][0], 1, xall_send_right, nbr_right[0], tag, &rho[0][0][0], 1, xall_recv_left, nbr_left[0], MPI_ANY_TAG, comm_cart, MPI_STATUS_IGNORE);

  //  MPI_Barrier();

  MPI_Sendrecv(&rho[0][0][0], 1, yall_send_left, nbr_left[1], tag, &rho[0][0][0], 1, yall_recv_right, nbr_right[1], MPI_ANY_TAG, comm_cart, MPI_STATUS_IGNORE);
  MPI_Sendrecv(&rho[0][0][0], 1, yall_send_right, nbr_right[1], tag, &rho[0][0][0], 1, yall_recv_left, nbr_left[1], MPI_ANY_TAG, comm_cart, MPI_STATUS_IGNORE);

  //  MPI_Barrier();

  MPI_Sendrecv(&rho[0][0][0], 1, zall_send_left, nbr_left[2], tag, &rho[0][0][0], 1, zall_recv_right, nbr_right[2], MPI_ANY_TAG, comm_cart, MPI_STATUS_IGNORE);
  MPI_Sendrecv(&rho[0][0][0], 1, zall_send_right, nbr_right[2], tag, &rho[0][0][0], 1, zall_recv_left, nbr_left[2], MPI_ANY_TAG, comm_cart, MPI_STATUS_IGNORE);

  MPI_Barrier(comm_cart);
}

///Gather NOT WORKING YET!!!!!!!!!!!!!!!!
void MHD::gather(){
  // int n;
  // MPI_Reduce(&NJ,&n,1,MPI_INT, MPI_SUM,0,comm_cart);
  // if(rank==0)
  //   cout<<"sum of NJ = "<<n/9<<endl;
  //
  int local_size_with_ghost[3] = {NITf,NJTf,NKTf};
  int local_sub_size[3] = {NI,NJ,NK};
  int local_start_idx[3]  = {ial,jal,kal};

  MPI_Datatype sub_type;
  MPI_Type_create_subarray(3, local_size_with_ghost, local_sub_size, local_start_idx, MPI_ORDER_C,MPI_FLOAT, &sub_type);
  MPI_Type_commit(&sub_type);

  int global_size[3] = {ni_global,nj_global,nk_global};
  int global_sub_size[3] = {NI,NJ,NK};
  int global_start_idx[3];
  global_start_idx[0] = NI*coord[0];
  global_start_idx[1] = NJ*coord[1];
  global_start_idx[2] = NK*coord[2];

  //  cout<<"rank = "<<rank<<" idx = "<<global_start_idx[0] <<" "<<global_start_idx[1]<<" "<<global_start_idx[2] <<endl;

  MPI_Datatype global_type;
  MPI_Type_create_subarray(3, global_size, global_sub_size, global_start_idx, MPI_ORDER_C,MPI_FLOAT, &global_type);
  MPI_Type_commit(&global_type);

  float p_g[ni_global][nj_global][nk_global];
  for(int i=0;i<ni_global;i++)
    for(int j=0;j<nj_global;j++)
      for(int k=0;k<nk_global;k++)
	p_g[i][j][k]=0.;


  int rcnt[size],disp[size];
  for(int i=0;i<size;i++){
    rcnt[i] = 1;
    disp[i] = 0;
    //    disp[i] = global_start_idx[0]*nj_global*nk_global + global_start_idx[1]*nk_global + global_start_idx[2];
  }
  //  MPI_Gatherv(&p[0][0][0], 1, sub_type, &p_g[0][0][0], rcnt, disp, global_type, 0, comm_cart);
  int tag=1;

  //  MPI_Sendrecv(&p[0][0][0],1,sub_type, rank, tag, &p_g[0][0][0],1,global_type, rank,MPI_ANY_TAG,comm_cart, MPI_STATUS_IGNORE);

  if(rank==0)
    MPI_Recv(&p_g[8][0][0],1,global_type, 18, MPI_ANY_TAG,comm_cart, MPI_STATUS_IGNORE);
  if(rank==18)
    MPI_Send(&p[0][0][0],1,sub_type, 0, tag, comm_cart);
  MPI_Barrier(comm_cart);

  if(rank==0)
    for(int i=0; i<ni_global; i++){
      for(int j=0; j<nj_global; j++){
	for(int k=0; k<nk_global; k++){
  	  cout<<p_g[i][j][k]<<" ";
	}
  	cout<<endl;
      }
      cout<<endl;
    }


  //  int n=1;
  //  int x[27];
  //  MPI_Gather(&n, 1, MPI_INT, x, 1, MPI_INT, 0, comm_cart);

  MPI_Type_free(&sub_type);
  MPI_Type_free(&global_type);

  MPI_Barrier(comm_cart);
}
