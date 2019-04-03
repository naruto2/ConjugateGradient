SRCS = main.c++ sparse.h crs.h operator.h solver.h \
	ConjugateGradient.h cgs.h bicgstab.h gmres.h bicg.h qmr.h

a.out: $(SRCS)
	c++ -Wall -DnoGPU main.c++ 

gpu: $(SRCS)
	ln -sf main.c++ main.cu
	nvcc main.cu -lcusparse -lcublas -o gpu

clean: ;
	rm -f ./a.out ./gpu ./main.cu
