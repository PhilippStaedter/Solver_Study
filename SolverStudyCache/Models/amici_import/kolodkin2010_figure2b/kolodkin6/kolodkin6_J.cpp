#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_kolodkin6(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[8] = 1.0*dwdx0;
    J[9] = 1.0*dwdx2 - 1.0*dwdx3;
    J[10] = -1.0*dwdx4;
    J[12] = -1.0*dwdx9;
    J[17] = -1.0*dwdx3;
    J[18] = -1.0*dwdx4 + 3.4444444444444402*dwdx5;
    J[20] = -1.0*dwdx9;
    J[21] = 3.4444444444444402*dwdx12;
    J[24] = 1.0*dwdx1;
    J[27] = 1.0*dwdx6 - 1.0*dwdx7;
    J[28] = -1.0*dwdx10;
    J[29] = 1.0*dwdx11;
    J[33] = 1.0*dwdx3;
    J[34] = 1.0*dwdx4;
    J[35] = 3.4444444444444402*dwdx7;
    J[36] = 3.4444444444444402*dwdx10 - 1.0*dwdx8 + 1.0*dwdx9;
    J[38] = -1.0*dwdx13;
    J[39] = -1.0*dwdx14;
    J[40] = -1.0*dwdx1;
    J[42] = -1.0*dwdx5;
    J[43] = -1.0*dwdx6;
    J[45] = -1.0*dwdx11 - 1.0*dwdx12;
    J[52] = -1.0*dwdx8;
    J[54] = -1.0*dwdx13;
    J[55] = -1.0*dwdx14;
    J[60] = 1.0*dwdx8;
    J[62] = 1.0*dwdx13;
    J[63] = 1.0*dwdx14;
}