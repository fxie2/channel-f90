F90 = gfortran
INC = /usr/include/

OBJ = typedef.o rbmat.o ffts.o dnsdata.o channel.o
#flags =  -malign-double -fall-intrinsics -ffree-line-length-none -I$(INC) -cpp  -g -fcheck=all -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -lfftw3
flags =  -malign-double -fall-intrinsics -ffree-line-length-none  -I$(INC) -cpp -lfftw3  -Ofast

channel: $(OBJ)
	$(F90) $(flags) -o $@ $(OBJ)
%.o : %.f90
	$(F90) $(flags) -c $<
clean: 
	rm *.mod *.o
