#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_hornberg1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[7] = -1.0*dwdx0;
    JB[8] = 1.0*dwdx0;
    JB[10] = 1.0*dwdx1;
    JB[11] = -1.0*dwdx1;
    JB[12] = 1.0*dwdx2;
    JB[13] = -1.0*dwdx2;
    JB[19] = -1.0*dwdx3;
    JB[20] = 1.0*dwdx3;
    JB[30] = 1.0*dwdx4;
    JB[31] = -1.0*dwdx4;
    JB[39] = -1.0*dwdx5;
    JB[40] = 1.0*dwdx5;
    JB[41] = 1.0*dwdx6;
    JB[42] = -1.0*dwdx6;
    JB[50] = 1.0*dwdx7;
    JB[51] = -1.0*dwdx7;
    JB[59] = -1.0*dwdx8;
    JB[60] = 1.0*dwdx8;
    JB[61] = 1.0*dwdx9;
    JB[62] = -1.0*dwdx9;
    JB[70] = 1.0*dwdx10;
    JB[71] = -1.0*dwdx10;
    JB[79] = -1.0*dwdx11;
    JB[80] = 1.0*dwdx11;
}