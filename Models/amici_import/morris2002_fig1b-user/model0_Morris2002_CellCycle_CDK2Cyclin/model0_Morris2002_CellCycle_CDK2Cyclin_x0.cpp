#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_Morris2002_CellCycle_CDK2Cyclin(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[2] = 9.9999999999999995e-8;
    x0[3] = 3.9999999999999998e-7;
}