#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_Singh2006(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 8.0;
    x0[7] = 100.0;
    x0[8] = 1000.0;
    x0[21] = 50.0;
    x0[28] = 60.0;
    x0[32] = 85.0;
    x0[33] = 19000.0;
    x0[36] = 67.0;
    x0[40] = 67.0;
    x0[41] = 41667.0;
    x0[45] = 67.0;
    x0[47] = 35000.0;
    x0[52] = 34.0;
    x0[57] = 0.80000000000000004;
    x0[59] = 12.0;
    x0[66] = 16667.0;
}