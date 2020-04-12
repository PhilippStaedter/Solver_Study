#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model5_panteleev3(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0;
    JB[5] = 1.0*dwdx0;
    JB[6] = -1.0*dwdx0;
    JB[9] = 1.0*dwdx1 - 1.0*dwdx2 + 1.0*dwdx3;
    JB[10] = -1.0*dwdx1;
    JB[11] = 1.0*dwdx2;
    JB[12] = 1.0*dwdx1;
    JB[13] = -1.0*dwdx2;
    JB[14] = 1.0*dwdx3;
    JB[15] = -1.0*dwdx3;
    JB[17] = 1.0*dwdx4;
    JB[18] = -1.0*dwdx4 + 1.0*dwdx5;
    JB[19] = -1.0*dwdx5;
    JB[20] = 1.0*dwdx4;
    JB[25] = -1.0*dwdx6;
    JB[27] = 1.0*dwdx6;
    JB[29] = -1.0*dwdx6;
    JB[33] = 1.0*dwdx7;
    JB[34] = -1.0*dwdx7;
    JB[36] = 1.0*dwdx7;
    JB[40] = 1.0*dwdx9;
    JB[41] = -1.0*dwdx8;
    JB[43] = 1.0*dwdx8;
    JB[45] = -1.0*dwdx8 + 1.0*dwdx9;
    JB[46] = -1.0*dwdx9;
    JB[48] = 1.0*dwdx10;
    JB[49] = 1.0*dwdx11;
    JB[53] = 1.0*dwdx10;
    JB[54] = -1.0*dwdx10 + 1.0*dwdx11;
    JB[55] = -1.0*dwdx11;
    JB[57] = 1.0*dwdx12;
    JB[62] = 1.0*dwdx12;
    JB[63] = -1.0*dwdx12;
}