#include "MHD.h"
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
#include <hdf5.h>
void MHD::hdfsave(const char* filename){
  hid_t err;
  hid_t dataspace, memspace, dataset;

  hid_t string_type;

  hid_t fap_id;

  /* define an info object to store MPI-IO information */
  MPI_Info finfo;

  char comment[] = "HDF5 file for MHD";

  hid_t file_id;

  int ierr;
  //  MPI_Status stat;

  /* create an MPI_INFO object -- on some platforms it is useful to
     pass some information onto the underlying MPI_File_open call */
  MPI_Info_create(&finfo);
  MPI_Info_set(finfo,(char*)"IBM_largeblock_io",(char*)"true");

  /* set the file access template for parallel IO access */
  fap_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fap_id, comm_cart, finfo);

  hid_t dist_id;
  dist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dist_id, H5FD_MPIO_COLLECTIVE);

  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fap_id);


  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, strlen(comment));
  int num_dims;
  hsize_t data_dimens_1d;
  num_dims = 1;
  data_dimens_1d = 1;
  dataspace = H5Screate_simple(num_dims, &data_dimens_1d, NULL);
  dataset = H5Dcreate(file_id, "comment", string_type, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  //  if (rank == 0) {
    err = H5Dwrite(dataset, string_type, H5S_ALL, dataspace, 
		   dist_id, comment);
    //  }
  H5Sclose(dataspace);
  H5Dclose(dataset);

  /*--------------------------------------------------------------------------
   * Now store the data array -- interior zones only
   *-------------------------------------------------------------------------*/
  num_dims = 3;
  // dataspace dimensions in total
  hsize_t data_dimens_3d[3];
  data_dimens_3d[0] = ni_global+1;
  data_dimens_3d[1] = nj_global+1;
  data_dimens_3d[2] = nk_global+1;

  /* figure out the offset into the dataspace for the current processor */    
  hsize_t data_start_3d[3];
  hsize_t stride_3d[3];
  hsize_t count_3d[3];
  data_start_3d[0] = offset_global[0];
  data_start_3d[1] = offset_global[1];
  data_start_3d[2] = offset_global[2];
  // stride for every proc
  stride_3d[0] = 1;
  stride_3d[1] = 1;
  stride_3d[2] = 1;
  // how many data in each dimension; for each proc, the counts are the same for memory and dataspace
  count_3d[0] = NI;
  count_3d[1] = NJ;
  count_3d[2] = NK;
  // store 1 more data if on the high boundary
  if(ihibound)
    count_3d[0] = NI+1;
  if(jhibound)
    count_3d[1] = NJ+1;
  if(khibound)
    count_3d[2] = NK+1;

  // data size on memory for each proc
  hsize_t mem_dimens_3d[3];
  mem_dimens_3d[0] = NITf;
  mem_dimens_3d[1] = NJTf;
  mem_dimens_3d[2] = NKTf;
  // start idx on memory for each proc
  hsize_t mem_start_3d[3];
  mem_start_3d[0] = NO2;
  mem_start_3d[1] = NO2;
  mem_start_3d[2] = NO2;

  // density ///
  /* create the dataspace (how it will be stored on disk).  Here we need
     to specify enough space for all of the processors' data.  Then we
     use the hyperslab to point to the portion of the whole dataspace
     that will store the current processor's data. */
  dataspace = H5Screate_simple(num_dims, data_dimens_3d, NULL);
  err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, 
			    data_start_3d, stride_3d, count_3d, NULL);
  /* now we need to create a memory space.  In this case, we want to
     extract the interior of the patch out when we write, so we'll need
     to use a hyperslab */
  memspace = H5Screate_simple(num_dims, mem_dimens_3d, NULL);
  //count_3d is the same as defined above
  /* now, use the HDF5 hyperslab function to pick out the portion of the data
     array that we intend to store */
  err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, 
                            mem_start_3d, stride_3d, count_3d, NULL);
  /* now create the dataset -- if it is the first time through, otherwise
     open the exisiting dataset */
  dataset = H5Dcreate(file_id, "x", H5T_NATIVE_FLOAT,
		      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* and finally write the data to disk -- in this call, we need to 
     include both the memory space and the data space. */
  err = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, 
		 dist_id, &(x[0][0][0]));
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  // pressure ///
  dataspace = H5Screate_simple(num_dims, data_dimens_3d, NULL);
  err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, 
  			    data_start_3d, stride_3d, count_3d, NULL);
  memspace = H5Screate_simple(num_dims, mem_dimens_3d, NULL);
  err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, 
                            mem_start_3d, stride_3d, count_3d, NULL);
  dataset = H5Dcreate(file_id, "y", H5T_NATIVE_FLOAT,
  		      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  err = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, 
  		 dist_id, &(y[0][0][0]));
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  //
  dataspace = H5Screate_simple(num_dims, data_dimens_3d, NULL);
  err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, data_start_3d, stride_3d, count_3d, NULL);

  memspace = H5Screate_simple(num_dims, mem_dimens_3d, NULL);
  err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start_3d, stride_3d, count_3d, NULL);

  dataset = H5Dcreate(file_id, "z", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  err = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, dist_id, &(z[0][0][0]));

  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  /*--------------------------------------------------------------------------
   * finally, close the file
   *-------------------------------------------------------------------------*/
  /* release the file access template */
  ierr = H5Pclose(fap_id);
  ierr = MPI_Info_free(&finfo);
  
  H5Fclose(file_id);

  MPI_Barrier(comm_cart);
}


void MHD::hdfread(const char* filename){
  hid_t err;
  hid_t dataspace, memspace, dataset;
  hid_t fap_id;
  hid_t file_id;
  MPI_Info finfo;

  int ierr;

  //
  fap_id = H5Pcreate(H5P_FILE_ACCESS);
  // ierr = H5Pset_sieve_buf_size(fap_id, 262144/2); 
  // ierr = H5Pset_alignment(fap_id, 524288/2, 262144/2);

  ierr = MPI_Info_create(&finfo);
  MPI_Info_set(finfo,(char*)"IBM_largeblock_io",(char*)"true");
  // ierr = MPI_Info_set(finfo, (char*)"access_style", (char*)"read_once");
  // ierr = MPI_Info_set(finfo, (char*)"collective_buffering", (char*)"true");
  // ierr = MPI_Info_set(finfo, (char*)"cb_block_size", (char*)"1048576");
  // ierr = MPI_Info_set(finfo, (char*)"cb_buffer_size", (char*)"4194304");

  ierr = H5Pset_fapl_mpio(fap_id, comm_cart, finfo);

  hid_t dist_id;
  dist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dist_id, H5FD_MPIO_COLLECTIVE);

  //
  file_id = H5Fopen(filename, H5F_ACC_RDONLY, fap_id);
  ierr = H5Pclose(fap_id);
  ierr = MPI_Info_free(&finfo);

  /////
  int num_dims;
  num_dims = 3;
  // dataspace dimensions in total
  hsize_t data_dimens_3d[3];
  data_dimens_3d[0] = ni_global+1;
  data_dimens_3d[1] = nj_global+1;
  data_dimens_3d[2] = nk_global+1;

  /* figure out the offset into the dataspace for the current processor */    
  hsize_t data_start_3d[3];
  hsize_t stride_3d[3];
  hsize_t count_3d[3];
  data_start_3d[0] = offset_global[0];
  data_start_3d[1] = offset_global[1];
  data_start_3d[2] = offset_global[2];
  stride_3d[0] = 1;
  stride_3d[1] = 1;
  stride_3d[2] = 1;
  // read 1 more for each
  count_3d[0] = NI+1;
  count_3d[1] = NJ+1;
  count_3d[2] = NK+1;
  // count_3d[0] = NI;
  // count_3d[1] = NJ;
  // count_3d[2] = NK;
  // if(ihibound)
  //   count_3d[0] = NI+1;
  // if(jhibound)
  //   count_3d[1] = NJ+1;
  // if(khibound)
  //   count_3d[2] = NK+1;

  // data size on memory for each proc
  hsize_t mem_dimens_3d[3];
  mem_dimens_3d[0] = NITf;
  mem_dimens_3d[1] = NJTf;
  mem_dimens_3d[2] = NKTf;
  // start idx on memory for each proc
  hsize_t mem_start_3d[3];
  mem_start_3d[0] = NO2;
  mem_start_3d[1] = NO2;
  mem_start_3d[2] = NO2;

  // density ///
  dataset = H5Dopen(file_id,"x",H5P_DEFAULT);
  //  dataset = H5Dopen(file_id,"/X_grid",H5P_DEFAULT);
  //  hid_t       datatype;
  //  datatype  = H5Dget_type(dataset);     /* datatype handle */ 
  dataspace = H5Dget_space(dataset);
  err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, 
			    data_start_3d, stride_3d, count_3d, NULL);
  memspace = H5Screate_simple(num_dims, mem_dimens_3d, NULL);
  //count_3d is the same as defined above
  err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, 
                            mem_start_3d, stride_3d, count_3d, NULL);
  err = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, 
		dist_id, &(x[0][0][0]));
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  // pressure ///
  dataset = H5Dopen(file_id,"y",H5P_DEFAULT);
  //  dataset = H5Dopen(file_id,"/Y_grid",H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);
  err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, 
			    data_start_3d, stride_3d, count_3d, NULL);
  memspace = H5Screate_simple(num_dims, mem_dimens_3d, NULL);
  err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, 
                            mem_start_3d, stride_3d, count_3d, NULL);
  err = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, 
		dist_id, &(y[0][0][0]));
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  //
  dataset = H5Dopen(file_id,"z",H5P_DEFAULT);
  //  dataset = H5Dopen(file_id,"/Z_grid",H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);
  err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, 
			    data_start_3d, stride_3d, count_3d, NULL);
  memspace = H5Screate_simple(num_dims, mem_dimens_3d, NULL);
  err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, 
                            mem_start_3d, stride_3d, count_3d, NULL);
  err = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, 
		dist_id, &(z[0][0][0]));
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);


  H5Fclose(file_id);

  MPI_Barrier(comm_cart);
}
void MHD::hdfdump(const char* filename){
  hid_t err;
  hid_t dataspace, memspace, dataset;

  hid_t string_type;

  hid_t fap_id;

  /* define an info object to store MPI-IO information */
  MPI_Info finfo;

  //  char *filename = "dump.hdf5";
  //  char filename[] = "mhd.hdf5";
  char comment[] = "HDF5 file for MHD";

  hid_t file_id;

  int ierr;
  //  MPI_Status stat;

  /*--------------------------------------------------------------------------
   * Now open the HDF5 file for output and write out the data array -- all
   * processors do the openning
   *-------------------------------------------------------------------------*/

  
  /* H5Fcreate takes several arguments in addition to the filename.  We 
     specify H5F_ACC_TRUNC to the second argument to tell it to overwrite 
     an existing file by the same name if it exists.  The next two 
     arguments are the file creation property list and the file access
     property lists.  These are used to pass options to the library about
     how to create the file, and how it will be accessed (ex. via mpi-io). */


  /* ---------------------------------------------------------------------
     platform dependent code goes here -- the access template must be
     tuned for a particular filesystem blocksize.  some of these 
     numbers are guesses / experiments, others come from the file system
     documentation.
     
     The sieve_buf_size should be equal a multiple of the disk block size
     ---------------------------------------------------------------------- */
  /* set the file access template for parallel IO access */
  fap_id = H5Pcreate(H5P_FILE_ACCESS);
  // ierr = H5Pset_sieve_buf_size(fap_id, 262144); 
  // ierr = H5Pset_alignment(fap_id, 524288, 262144);

  /* create an MPI_INFO object -- on some platforms it is useful to
     pass some information onto the underlying MPI_File_open call */

  ierr = MPI_Info_create(&finfo);
  // //in C++, string literal "..." is const char*, but the function need char*
  MPI_Info_set(finfo,(char*)"IBM_largeblock_io",(char*)"true");
  // ierr = MPI_Info_set(finfo, (char*)"access_style", (char*)"write_once");
  // ierr = MPI_Info_set(finfo, (char*)"collective_buffering", (char*)"true");
  // ierr = MPI_Info_set(finfo, (char*)"cb_block_size", (char*)"1048576");
  // ierr = MPI_Info_set(finfo, (char*)"cb_buffer_size", (char*)"4194304");
  /* ----------------------------------------------------------------------
     end of platform dependent properties
     ---------------------------------------------------------------------- */

  /* tell the HDF5 library that we want to use MPI-IO to do the writing */
  ierr = H5Pset_fapl_mpio(fap_id, comm_cart, finfo);

  hid_t dist_id;
  dist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dist_id, H5FD_MPIO_COLLECTIVE);

  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, 
			      fap_id);

  /* release the file access template */
  ierr = H5Pclose(fap_id);
  ierr = MPI_Info_free(&finfo);

  /*--------------------------------------------------------------------------
   * Start by creating a record that stores the comment -- this will only
   * be written by processor 0, but all procs create the dataspace.
   *-------------------------------------------------------------------------*/
  
  /* to store a string in HDF5, we need to build a special type out of 
     character types.  We use H5T_C_S1 -- a one-byte, null-terminated
     string -- as the basis for our type */
  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, strlen(comment));

  /* next we create a dataspace -- these describes how the data is stored
     in the file.  We need to tell it the number of dims and the dimensions of the
     record we are storing */
  int num_dims;
  hsize_t data_dimens_1d;
  num_dims = 1;
  data_dimens_1d = 1;
  /* The data we are storing will be contiguous on disk, so we can just 
     create a simple dataspace.  H5Screate_simple takes 3 arguments, the
     rank, an array (or size rank) of the dimensions, and a maxdims argument
     that can be used to create unlimited sized datasets -- we use NULL 
     here */
  dataspace = H5Screate_simple(num_dims, &data_dimens_1d, NULL);
  
  /* next we create the dataset -- this tells the library the type of 
     data we are storing, a record name that we will refer to it by,
     and the dataspace that describes the data.  The last argument 
     refers to the dataset property list -- this can be used to tell
     the library to use compression, or various properties about how
     it is stored on disk.  We accept the defaults here */
  /* This works on babylon5:*/
  dataset = H5Dcreate(file_id, "comment", string_type, 
		      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* finally we write it to the file.  We give H5Dwrite the data set to write
     to, the type of data, information on how the data is stored in memory
     (in this case it is trivial, so we use H5S_ALL), the dataspace, a
     data transfer property list (we take the default), and a pointer to the 
     data buffer */
  //  if (rank == 0) {
    err = H5Dwrite(dataset, string_type, H5S_ALL, dataspace, 
		   dist_id, comment);
    //  }

  /* after the write, we are done with this dataspace and dataset, so we
     need to free the space we allocated for them. */
  H5Sclose(dataspace);
  H5Dclose(dataset);

  /*--------------------------------------------------------------------------
   * Now store the data array -- interior zones only
   *-------------------------------------------------------------------------*/
  num_dims = 3;
  // dataspace dimensions in total
  hsize_t data_dimens_3d[3];
  data_dimens_3d[0] = ni_global+1;
  data_dimens_3d[1] = nj_global+1;
  data_dimens_3d[2] = nk_global+1;

  /* figure out the offset into the dataspace for the current processor */    
  hsize_t data_start_3d[3];
  hsize_t stride_3d[3];
  hsize_t count_3d[3];
  data_start_3d[0] = offset_global[0];
  data_start_3d[1] = offset_global[1];
  data_start_3d[2] = offset_global[2];
  // stride for every proc
  stride_3d[0] = 1;
  stride_3d[1] = 1;
  stride_3d[2] = 1;
  // how many data in each dimension; for each proc, the counts are the same for memory and dataspace
  count_3d[0] = NI;
  count_3d[1] = NJ;
  count_3d[2] = NK;
  // store 1 more data if on the high boundary
  if(ihibound)
    count_3d[0] = NI+1;
  if(jhibound)
    count_3d[1] = NJ+1;
  if(khibound)
    count_3d[2] = NK+1;

  // data size on memory for each proc
  hsize_t mem_dimens_3d[3];
  mem_dimens_3d[0] = NITf;
  mem_dimens_3d[1] = NJTf;
  mem_dimens_3d[2] = NKTf;
  // start idx on memory for each proc
  hsize_t mem_start_3d[3];
  mem_start_3d[0] = NO2;
  mem_start_3d[1] = NO2;
  mem_start_3d[2] = NO2;

  // density ///
  /* create the dataspace (how it will be stored on disk).  Here we need
     to specify enough space for all of the processors' data.  Then we
     use the hyperslab to point to the portion of the whole dataspace
     that will store the current processor's data. */
  dataspace = H5Screate_simple(num_dims, data_dimens_3d, NULL);
  err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, 
			    data_start_3d, stride_3d, count_3d, NULL);
  /* now we need to create a memory space.  In this case, we want to
     extract the interior of the patch out when we write, so we'll need
     to use a hyperslab */
  memspace = H5Screate_simple(num_dims, mem_dimens_3d, NULL);
  //count_3d is the same as defined above
  /* now, use the HDF5 hyperslab function to pick out the portion of the data
     array that we intend to store */
  err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, 
                            mem_start_3d, stride_3d, count_3d, NULL);
  /* now create the dataset -- if it is the first time through, otherwise
     open the exisiting dataset */
  dataset = H5Dcreate(file_id, "x", H5T_NATIVE_FLOAT,
		      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* and finally write the data to disk -- in this call, we need to 
     include both the memory space and the data space. */
  err = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, 
		 dist_id, &(x[0][0][0]));
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  // pressure ///
  dataspace = H5Screate_simple(num_dims, data_dimens_3d, NULL);
  err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, 
  			    data_start_3d, stride_3d, count_3d, NULL);
  memspace = H5Screate_simple(num_dims, mem_dimens_3d, NULL);
  err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, 
                            mem_start_3d, stride_3d, count_3d, NULL);
  dataset = H5Dcreate(file_id, "y", H5T_NATIVE_FLOAT,
  		      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  err = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, 
  		 dist_id, &(y[0][0][0]));
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  //
  dataspace = H5Screate_simple(num_dims, data_dimens_3d, NULL);
  err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, 
  			    data_start_3d, stride_3d, count_3d, NULL);
  memspace = H5Screate_simple(num_dims, mem_dimens_3d, NULL);
  err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, 
                            mem_start_3d, stride_3d, count_3d, NULL);
  dataset = H5Dcreate(file_id, "z", H5T_NATIVE_FLOAT,
  		      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  err = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, 
  		 dist_id, &(z[0][0][0]));
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  /*--------------------------------------------------------------------------
   * finally, close the file
   *-------------------------------------------------------------------------*/
  
  H5Fclose(file_id);

  MPI_Barrier(comm_cart);
}
