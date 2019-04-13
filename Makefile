SRCS = main.c++ solver.h 

a.out: $(SRCS)
	c++ -Wall -DnoGPU -I. main.c++ src/mmio.c++ 

gpu: $(SRCS)
	ln -sf main.c++ main.cu
	c++ -c -I. src/mmio.c++
	nvcc main.cu mmio.o -lcusparse -lcublas -o gpu

clean: ;
	rm -f ./a.out ./gpu ./main.cu ./yosen.o ./honsen.o ./mmio.o

speedtest: gpu
	export PSC98=0; time ./gpu
	export PSC98=1; time ./gpu
	export PSC98=2; time ./gpu
	export PSC98=3; time ./gpu
	export PSC98=4; time ./gpu
	export PSC98=5; time ./gpu

