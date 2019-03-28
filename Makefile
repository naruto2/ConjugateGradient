a.out: main.c++ ConjugateGradient.h cgs.h crs.h sparse.h operator.h
	c++ -DnoGPU main.c++ 

gpu: cgs.h
	ln -sf main.c++ main.cu
	nvcc main.cu -lcusparse -lcublas -o gpu

clean: ;
	rm -f ./a.out ./gpu ./main.cu
