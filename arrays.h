#include <cstring>
void Delete1DArray(void *m);
void Delete2DArray(void **m);
void Delete3DArray(void ***m);
void Delete4DArray(void ****m);

char * New1DArray(int ni,size_t dszie);
char **New2DArray(int ni, int nj, size_t dsize);
char ***New3DArray(int ni, int nj, int nk, size_t dsize);
char ****New4DArray(int ns, int ni, int nj, int nk, size_t dsize);
