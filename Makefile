a.out: main_sample.cu ConjugateGradient.h csrmatrix.h sparsematrix.h
	nvcc main_sample.cu -lcusparse -lcublas
