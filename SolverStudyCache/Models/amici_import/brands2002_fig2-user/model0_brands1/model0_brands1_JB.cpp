#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_brands1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0;
    JB[8] = -1.0*dwdx0;
    JB[22] = -1.0*dwdx2;
    JB[23] = -1.0*dwdx1;
    JB[24] = 1.0*dwdx1 + 1.0*dwdx2;
    JB[32] = -1.0*dwdx1;
    JB[66] = -1.0*dwdx3;
    JB[69] = -1.0*dwdx5;
    JB[71] = -1.0*dwdx5;
    JB[72] = 1.0*dwdx3 + 1.0*dwdx4 + 1.0*dwdx5 + 1.0*dwdx6;
    JB[73] = -1.0*dwdx4;
    JB[75] = -2.0*dwdx6;
    JB[76] = 1.0*dwdx3;
    JB[79] = -1.0*dwdx9;
    JB[80] = -1.0*dwdx8;
    JB[82] = -1.0*dwdx8;
    JB[83] = -1.0*dwdx7;
    JB[84] = 1.0*dwdx7 + 1.0*dwdx8 + 1.0*dwdx9;
    JB[87] = 1.0*dwdx9;
    JB[100] = -1.0*dwdx10;
    JB[103] = -1.0*dwdx10;
    JB[108] = 1.0*dwdx10;
    JB[110] = -1.0*dwdx11;
    JB[112] = -1.0*dwdx12;
    JB[116] = 1.0*dwdx11;
    JB[117] = 1.0*dwdx12;
    JB[120] = 1.0*dwdx11 + 1.0*dwdx12;
}