FC=gfortran
FLAGS= -g -O0 -ffpe-trap=invalid,zero,overflow

FOBJ = types.o const.o utils.o particles.o output.o kepler.o drift.o  wh91.o kick.o step.o aside.o

%.o: %.F90
	${FC} -c ${FLAGS} $<

all: ${FOBJ}
	${FC} -o aside.x ${FLAGS} ${FOBJ}

clean:
	rm *.mod *.o