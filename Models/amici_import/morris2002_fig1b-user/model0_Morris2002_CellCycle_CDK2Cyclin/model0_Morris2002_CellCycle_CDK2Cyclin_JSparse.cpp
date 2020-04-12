#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_model0_Morris2002_CellCycle_CDK2Cyclin(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = -1000000000000.0*dwdx0 + 1000000000000.0*dwdx1;
    JSparse[1] = 1000000000000.0*dwdx0;
    JSparse[2] = -1000000000000.0*dwdx1;
    JSparse[3] = -1000000000000.0*dwdx1;
    JSparse[4] = -1000000000000.0*dwdx2;
    JSparse[5] = 1000000000000.0*dwdx2;
    JSparse[6] = 1000000000000.0*dwdx3;
    JSparse[7] = -1000000000000.0*dwdx3;
    JSparse[8] = -1000000000000.0*dwdx3;
    JSparse[9] = 1000000000000.0*dwdx4;
    JSparse[10] = -1000000000000.0*dwdx4;
    JSparse[11] = -1000000000000.0*dwdx4;
}