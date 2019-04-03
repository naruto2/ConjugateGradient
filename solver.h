#ifndef SOLVER_H
#define SOLVER_H

#ifndef noGPU
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>
#include <thrust/device_vector.h>
using namespace thrust;
#endif

#define SOLVERROR_NONE      0
#define SOLVERROR_MAXIT     1
#define SOLVERROR_BREAKDOWN 2

#include <iostream>
#include <cmath>
#include "sparse.h"
using namespace std;
using namespace sparse;
#include "crs.h"
#include "operator.h"
#include "ConjugateGradient.h"
#include "cgs.h"
#include "bicgstab.h"
#include "gmres.h"

#endif
