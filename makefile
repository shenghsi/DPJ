#ccc = g++
MPICXX = mpicxx
HDF5path = /usr

# headers to get .o files:
# By default, gcc searches the following directories for header files:
# /usr/local/include/
# /usr/include/
# use -I${HDF5path}/include to include headers that are not in the default directories.
# libraries for linking:
# and the following directories for libraries:
# /usr/local/lib/
# /usr/lib/
# use -L${HDF5path}/lib to specify the hdf path if the library is not in the default directories (/usr). otherwise just use -lhdf5

INC = -I${HDF5path}/include

# H5Dcreate may use different lib on different thayer computers
LIB = -L${HDF5path}/lib -lhdf5

flags = ${INC} -O2

OFILES= main.o MHD.o arrays.o global.o IO.o

all:	$(OFILES)  #declare the files needed if not exist yet
	$(MPICXX) $(OFILES) $(LIB) -o test

.C.o:
	$(MPICXX)  $(flags) -c $< -o $@  
#	$(MPICXX) -c $<	
clean:
	rm *.o test *~ *.hdf5
