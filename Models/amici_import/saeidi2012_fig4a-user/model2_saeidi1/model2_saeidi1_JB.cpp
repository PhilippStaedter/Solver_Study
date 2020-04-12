#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model2_saeidi1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0;
    JB[4] = -1.0*dwdx0;
    JB[11] = 1.0*dwdx1;
    JB[12] = -1.0*dwdx1;
    JB[17] = 1.0*dwdx1;
    JB[21] = 1.0*dwdx3;
    JB[22] = 1.0*dwdx2 - 1.0*dwdx3;
    JB[27] = 1.0*dwdx3;
    JB[28] = -1.0*dwdx2;
    JB[33] = 1.0*dwdx4 + 1.0*dwdx5;
    JB[35] = -1.0*dwdx4;
    JB[36] = 1.0*dwdx5;
    JB[37] = -1.0*dwdx5;
    JB[40] = 1.0*dwdx6;
    JB[43] = -1.0*dwdx7;
    JB[44] = -1.0*dwdx6 + 1.0*dwdx7;
    JB[63] = 1.0*dwdx8;
    JB[66] = 1.0*dwdx8 + 1.0*dwdx9;
    JB[67] = -1.0*dwdx8;
    JB[69] = -1.0*dwdx9;
    JB[71] = 1.0*dwdx11;
    JB[72] = -1.0*dwdx11;
    JB[73] = 1.0*dwdx10;
    JB[76] = 1.0*dwdx10;
    JB[77] = -1.0*dwdx10 + 1.0*dwdx11;
    JB[82] = 1.0*dwdx12;
    JB[88] = -1.0*dwdx12;
}