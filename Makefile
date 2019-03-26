a.out: main_sample.c++ ConjugateGradient.h crs.h simple.h custo.h printmatrix.h operator.h
	c++ -DnoGPU main_sample.c++ 

clean: ;
	rm -f ./a.out
