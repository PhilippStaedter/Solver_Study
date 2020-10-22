#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_Morris2002_CellCycle_CDK2Cyclin(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1000000000000.0*dwdx0 - 1000000000000.0*dwdx1;
    JB[1] = -1000000000000.0*dwdx0;
    JB[2] = 1000000000000.0*dwdx1;
    JB[3] = 1000000000000.0*dwdx1;
    JB[4] = 1000000000000.0*dwdx2;
    JB[5] = -1000000000000.0*dwdx2;
    JB[8] = -1000000000000.0*dwdx3;
    JB[10] = 1000000000000.0*dwdx3;
    JB[11] = 1000000000000.0*dwdx3;
    JB[12] = -1000000000000.0*dwdx4;
    JB[14] = 1000000000000.0*dwdx4;
    JB[15] = 1000000000000.0*dwdx4;
}