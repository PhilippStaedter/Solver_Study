#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_fung1_Fig3A_Vgly_0_5(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0 + 1.0*dwdx1;
    JB[1] = -1.0*dwdx0;
    JB[9] = 1.0*dwdx4;
    JB[10] = -1.0*dwdx2;
    JB[13] = -1.0*dwdx3;
    JB[14] = -1.0*dwdx4;
    JB[16] = -1.0*dwdx6;
    JB[18] = 1.0*dwdx5;
    JB[22] = 1.0*dwdx6;
    JB[27] = -1.0*dwdx7 + 1.0*dwdx8;
    JB[30] = 1.0*dwdx7;
    JB[35] = 1.0*dwdx9;
    JB[45] = 1.0*dwdx11;
    JB[47] = -1.0*dwdx10;
    JB[48] = -1.0*dwdx14;
    JB[49] = 1.0*dwdx13;
    JB[51] = -1.0*dwdx12;
    JB[54] = 1.0*dwdx12 - 1.0*dwdx13 + 1.0*dwdx14;
    JB[56] = 1.0*dwdx16;
    JB[57] = -1.0*dwdx16;
    JB[63] = 1.0*dwdx15;
}