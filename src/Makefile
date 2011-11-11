FC=gfortran
FLAGS= -g -O0 -ffpe-trap=invalid,zero,overflow -std=f95

FOBJ = types.o const.o utils.o particles.o input.o output.o kepler.o drift.o  wh91.o kick.o step.o corrector.o aside.o

%.o: %.f95
	${FC} -c ${FLAGS} $<

all: ${FOBJ}
	${FC} -o ../exe/aside.x ${FLAGS} ${FOBJ}

clean:
	rm *.mod *.o