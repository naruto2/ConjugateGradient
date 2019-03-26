a.out: main_sample.cu ConjugateGradient.h crs.h simple.h custo.h printmatrix.h operator.h
	nvcc main_sample.cu -lcusparse -lcublas

clean: ;
	rm -f ./a.out
