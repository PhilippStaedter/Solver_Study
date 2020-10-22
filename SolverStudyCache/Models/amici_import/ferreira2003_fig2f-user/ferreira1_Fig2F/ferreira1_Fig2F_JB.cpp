#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_ferreira1_Fig2F(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0 + 1.0*dwdx1;
    JB[1] = -1.0*dwdx1;
    JB[5] = -1.0*dwdx0;
    JB[14] = 1.0*dwdx2 + 1.0*dwdx3;
    JB[15] = -1.0*dwdx3;
    JB[16] = 1.0*dwdx2;
    JB[17] = -1.0*dwdx2;
    JB[19] = -1.0*dwdx4;
    JB[21] = 1.0*dwdx4 + 1.0*dwdx5;
    JB[22] = 1.0*dwdx4;
    JB[25] = -1.0*dwdx7;
    JB[26] = 1.0*dwdx6;
    JB[27] = 1.0*dwdx7;
    JB[28] = 1.0*dwdx6 + 1.0*dwdx7;
    JB[29] = -1.0*dwdx6;
    JB[30] = -1.0*dwdx9;
    JB[31] = -1.0*dwdx10;
    JB[32] = -1.0*dwdx8;
    JB[33] = -1.0*dwdx12;
    JB[34] = -1.0*dwdx11 - 1.0*dwdx8;
    JB[35] = 1.0*dwdx10 + 1.0*dwdx13 + 1.0*dwdx8 + 1.0*dwdx9;
}