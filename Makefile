SRCS = main.c++ sparse.h crs.h operator.h solver.h \
	ConjugateGradient.h cgs.h bicgstab.h gmres.h bicg.h qmr.h

a.out: $(SRCS)
	cc -c yosen.c
	c++ -Wall -DnoGPU main.c++ yosen.o

gpu: $(SRCS)
	ln -sf main.c++ main.cu
	cc -c yosen.c
	nvcc main.cu yosen.o -lcusparse -lcublas -o gpu

clean: ;
	rm -f ./a.out ./gpu ./main.cu ./yosen.o
