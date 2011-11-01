FC=gfortran
FLAGS= -g -O0 -ffpe-trap=invalid,zero,overflow

FOBJ = types.o const.o utils.o particles.o kepler.o wh91.o drift.o kick.o aside.o

%.o: %.F90
	${FC} -c ${FLAGS} $<

all: ${FOBJ}
	${FC} -o aside.x ${FLAGS} ${FOBJ}

clean:
	rm *.mod *.o