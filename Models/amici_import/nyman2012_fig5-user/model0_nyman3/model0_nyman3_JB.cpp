#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_nyman3(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0;
    JB[5] = -1.0*dwdx0;
    JB[10] = 1.0*dwdx1;
    JB[11] = -1.0*dwdx1;
    JB[19] = -1.0*dwdx3;
    JB[20] = 1.0*dwdx3;
    JB[25] = 1.0*dwdx2;
    JB[26] = -1.0*dwdx2;
    JB[27] = -1.0*dwdx4;
    JB[30] = 1.0*dwdx4;
    JB[37] = 1.0*dwdx6;
    JB[38] = -1.0*dwdx6;
    JB[39] = -1.0*dwdx5;
    JB[40] = 1.0*dwdx5;
    JB[45] = -1.0*dwdx7;
    JB[50] = 1.0*dwdx7 + 1.0*dwdx8;
    JB[51] = -1.0*dwdx8;
    JB[54] = -1.0*dwdx10;
    JB[55] = 1.0*dwdx11;
    JB[56] = -1.0*dwdx11;
    JB[58] = -1.0*dwdx9;
    JB[60] = 1.0*dwdx10 + 1.0*dwdx9;
    JB[70] = 1.0*dwdx12;
    JB[71] = -1.0*dwdx12;
    JB[73] = 1.0*dwdx14;
    JB[74] = -1.0*dwdx14;
    JB[75] = -1.0*dwdx13;
    JB[76] = 1.0*dwdx13;
    JB[79] = -1.0*dwdx15;
    JB[80] = 1.0*dwdx15;
}