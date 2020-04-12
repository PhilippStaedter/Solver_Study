#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_model0_Morris2002_CellCycle_CDK2Cyclin(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = -1000000000000.0*flux_r0 + 1000000000000.0*flux_r1;
    xdot[1] = 1000000000000.0*flux_r0;
    xdot[2] = -1000000000000.0*flux_r1;
    xdot[3] = -1000000000000.0*flux_r1;
}