OCTFILES=basisfun.oct bspeval.oct nrb_srf_basisfun__.oct surfderivcpts.oct basisfunder.oct \
curvederivcpts.oct nrb_srf_basisfun_der__.oct surfderiveval.oct bspderiv.oct \
nrbsurfderiveval.oct tbasisfun.oct

all: $(OCTFILES)

low_level_functions.o: low_level_functions.cc
	mkoctfile -c $<

%.oct:  %.cc low_level_functions.o
	mkoctfile $< low_level_functions.o

clean:
	-rm -f *.o core octave-core *.oct *~
# DO NOT DELETE
