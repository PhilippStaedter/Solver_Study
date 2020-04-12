#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_Morris2002_CellCycle_CDK2Cyclin(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1000000000000.0*dwdx0 + 1000000000000.0*dwdx1;
    J[1] = -1000000000000.0*dwdx2;
    J[2] = 1000000000000.0*dwdx3;
    J[3] = 1000000000000.0*dwdx4;
    J[4] = 1000000000000.0*dwdx0;
    J[5] = 1000000000000.0*dwdx2;
    J[8] = -1000000000000.0*dwdx1;
    J[10] = -1000000000000.0*dwdx3;
    J[11] = -1000000000000.0*dwdx4;
    J[12] = -1000000000000.0*dwdx1;
    J[14] = -1000000000000.0*dwdx3;
    J[15] = -1000000000000.0*dwdx4;
}