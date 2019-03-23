a.out: main_sample.cu ConjugateGradient.h csrmatrix.h simple.h custo.h printmatrix.h
	nvcc main_sample.cu -lcusparse -lcublas

clean: ;
	rm -f ./a.out
