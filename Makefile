SRCS = main.c++ sparse.h crs.h operator.h solver.h \
	ConjugateGradient.h cgs.h bicgstab.h gmres.h bicg.h qmr.h

a.out: $(SRCS)
	cc -c yosen.c
	c++ -Wall -DnoGPU main.c++ mmio.c++ yosen.o

gpu: $(SRCS)
	ln -sf main.c++ main.cu
	cc -c yosen.c
	c++ -c mmio.c++
	nvcc main.cu mmio.o yosen.o -lcusparse -lcublas -o gpu

clean: ;
	rm -f ./a.out ./gpu ./main.cu ./yosen.o ./honsen.o ./mmio.o

speedtest: gpu
	export PSC98=0; time ./gpu
	export PSC98=1; time ./gpu
	export PSC98=2; time ./gpu
	export PSC98=3; time ./gpu
	export PSC98=4; time ./gpu
	export PSC98=5; time ./gpu

