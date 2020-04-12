#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_kholodenko2(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0;
    JB[1] = -1.0*dwdx0;
    JB[8] = -1.0*dwdx2;
    JB[9] = 1.0*dwdx1 + 1.0*dwdx2;
    JB[10] = -1.0*dwdx1;
    JB[17] = -1.0*dwdx4;
    JB[18] = 1.0*dwdx4;
    JB[20] = 1.0*dwdx3;
    JB[21] = -1.0*dwdx3;
    JB[27] = 1.0*dwdx5;
    JB[30] = -1.0*dwdx5;
    JB[36] = 1.0*dwdx6;
    JB[37] = -1.0*dwdx6;
    JB[43] = 1.0*dwdx8;
    JB[44] = -1.0*dwdx7;
    JB[45] = 1.0*dwdx7;
    JB[46] = -1.0*dwdx8 + 1.0*dwdx9;
    JB[47] = -1.0*dwdx9;
    JB[51] = -1.0*dwdx11;
    JB[54] = 1.0*dwdx10 + 1.0*dwdx11;
    JB[55] = -1.0*dwdx10;
    JB[56] = 1.0*dwdx13;
    JB[57] = -1.0*dwdx13 + 1.0*dwdx14;
    JB[58] = -1.0*dwdx14;
    JB[62] = -1.0*dwdx12;
    JB[63] = 1.0*dwdx12;
}