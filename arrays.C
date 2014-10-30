#include "arrays.h"
#include <stdlib.h>
#include <iostream>
using namespace std;

void print(const char *msg){
  cout<<msg<<endl;
}
void Delete1DArray(void *m){
  delete [] (char*)m;
}
void Delete2DArray(void **m){
  delete [] (char*)m[0];
  delete [] (char*)m;
}
void Delete3DArray(void ***m){
  delete [] (char*)m[0][0];
  delete [] (char*)m[0];
  delete [] (char*)m;
}
void Delete4DArray(void ****m){
  delete [] (char*)m[0][0][0];
  delete [] (char*)m[0][0];
  delete [] (char*)m[0];
  delete [] (char*)m;
}

char *New1DArray(int ni, size_t dsize){
  char *m = new char[ni*dsize];
  if(!m){
    cout<<"! Allocation failure in New1DArray"<<endl;
    exit(1);
  }
  return m;
}

char **New2DArray(int ni, int nj, size_t dsize){
  int i,j;
  char **m = new char*[ni];
  if(!m){
    cout<<"! Allocation failure in New2DArray"<<endl;
    exit(1);
  }

  m[0] = new char[ni*nj*dsize];
  if(!m[0]){
    cout<<"! Allocation failure in New2DArray"<<endl;
    exit(1);
  }

  for(i=1;i<ni;i++)
    m[i] = m[i-1] + nj*dsize; // for sizeof(type) = dsize, the address offset for the data is nj*dsize

  return m;
}

char ***New3DArray(int ni, int nj, int nk, size_t dsize){
  int i,j,k;
  char ***m = new char**[ni];
  if(!m){
    cout<<"! Allocation failure in New3DArray"<<endl;
    exit(1);
  }

  m[0] = new char*[ni*nj];
  if(!m[0]){
    cout<<"! Allocation failure in New3DArray"<<endl;
    exit(1);
  }

  m[0][0] = new char[ni*nj*nk*dsize];
  if(!m[0][0]){
    cout<<"! Allocation failure in New3DArray"<<endl;
    exit(1);
  }

  for(i=1;i<ni;i++)
    m[i] = m[i-1] + nj; // address offset for the address

  for(j=1;j<nj;j++)
    m[0][j] = m[0][j-1] + nk*dsize;

  for(i=1;i<ni;i++)
    m[i][0] = m[i-1][0] + nj*nk*dsize;

  for(i=1;i<ni;i++)
    for(j=1;j<nj;j++)
      m[i][j] = m[i][j-1] + nk*dsize;

  for(i=0;i<ni;i++)
    for(j=1;j<nj;j++){
      if(m[i][j]==NULL){
	cout<<"! Allocation failure in New3DArray"<<endl;
	exit(1);
      }
    }
  return m;
}

char ****New4DArray(int ni, int nj, int nk, int nl, size_t dsize){
  int i,j,k,l;
  char ****m = new char***[ni];
  if(!m){
    cout<<"! Allocation failure in New4DArray"<<endl;
    exit(1);
  }
  m[0] = new char**[ni*nj];
  if(!m[0]){
    cout<<"! Allocation failure in New4DArray"<<endl;
    exit(1);
  }
  m[0][0] = new char*[ni*nj*nk];
  if(!m[0][0]){
    cout<<"! Allocation failure in New4DArray"<<endl;
    exit(1);
  }
  m[0][0][0] = new char[ni*nj*nk*nl*dsize];
  if(!m[0][0][0]){
    cout<<"! Allocation failure in New4DArray"<<endl;
    exit(1);
  }

  for(i=1;i<ni;i++)
    m[i] = m[i-1] + nj; // address offset for the address
  //
  for(i=1;i<ni;i++)
    m[i][0] = m[i-1][0] + nj*nk;

  for(j=1;j<nj;j++)
    m[0][j] = m[0][j-1] + nk;

  for(i=1;i<ni;i++)
    for(j=1;j<nj;j++)
      m[i][j] = m[i][j-1] + nk;

/* ---------------------------
     triple subscript:
     
     (i,0,0) (0,j,0) (0,0,k)
     (i,j,0) (i,0,k) (0,j,k)
     (i,j,k)
   --------------------------- */

  for(i = 1; i < ni; i++) {
    m[i][0][0] = m[i-1][0][0] + nj*nk*nl*dsize;
  }
  for (j = 1; j < nj; j++) {
    m[0][j][0] = m[0][j - 1][0] + nk*nl*dsize;
  }
  for (k = 1; k < nk; k++) {
    m[0][0][k] = m[0][0][k - 1] + nl*dsize;
  }
  for (i = 1; i < ni; i++) {
    for (j = 1; j < nj; j++) {
      m[i][j][0] = m[i][j - 1][0] + nk*nl*dsize;
    }
  }
  for (i = 1; i < ni; i++) {
    for (k = 1; k < nk; k++) {
      m[i][0][k] = m[i][0][k - 1] + nl*dsize;
    }
  }
  for (j = 1; j < nj; j++) {
    for (k = 1; k < nk; k++) {
      m[0][j][k] = m[0][j][k - 1] + nl*dsize;
    }
  }
  for (i = 1; i < ni; i++) {
    for (j = 1; j < nj; j++) {
      for (k = 1; k < nk; k++) {
        m[i][j][k] = m[i][j][k - 1] + nl*dsize;
      }
    }
  }

 return m;
}
