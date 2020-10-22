#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_model0_Morris2002_CellCycle_CDK2Cyclin(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = 1000000000000.0*dwdx0 - 1000000000000.0*dwdx1;
    JSparseB[1] = 1000000000000.0*dwdx2;
    JSparseB[2] = -1000000000000.0*dwdx3;
    JSparseB[3] = -1000000000000.0*dwdx4;
    JSparseB[4] = -1000000000000.0*dwdx0;
    JSparseB[5] = -1000000000000.0*dwdx2;
    JSparseB[6] = 1000000000000.0*dwdx1;
    JSparseB[7] = 1000000000000.0*dwdx3;
    JSparseB[8] = 1000000000000.0*dwdx4;
    JSparseB[9] = 1000000000000.0*dwdx1;
    JSparseB[10] = 1000000000000.0*dwdx3;
    JSparseB[11] = 1000000000000.0*dwdx4;
}