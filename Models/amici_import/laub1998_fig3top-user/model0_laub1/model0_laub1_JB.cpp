#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_laub1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0;
    JB[5] = -1.0*dwdx1;
    JB[8] = 1.0*dwdx3;
    JB[9] = -1.0*dwdx2;
    JB[16] = 1.0*dwdx5;
    JB[19] = 1.0*dwdx4;
    JB[20] = 1.0*dwdx6;
    JB[22] = 1.0*dwdx8;
    JB[24] = 1.0*dwdx7;
    JB[28] = -1.0*dwdx10;
    JB[29] = -1.0*dwdx11;
    JB[32] = 1.0*dwdx9;
    JB[40] = 1.0*dwdx12;
    JB[41] = -1.0*dwdx13;
    JB[45] = 1.0*dwdx16;
    JB[46] = -1.0*dwdx14;
    JB[48] = 1.0*dwdx15;
}