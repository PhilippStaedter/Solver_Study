#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model2_panteleev3(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0;
    J[5] = -1.0*dwdx9;
    J[6] = -1.0*dwdx10;
    J[9] = -1.0*dwdx1 + 1.0*dwdx2 - 1.0*dwdx3;
    J[10] = -1.0*dwdx4;
    J[11] = 1.0*dwdx6;
    J[12] = -1.0*dwdx7;
    J[13] = 1.0*dwdx8;
    J[14] = -1.0*dwdx11;
    J[15] = -1.0*dwdx12;
    J[17] = 1.0*dwdx1;
    J[18] = 1.0*dwdx4 - 1.0*dwdx5;
    J[20] = 1.0*dwdx7;
    J[25] = -1.0*dwdx2;
    J[26] = 1.0*dwdx5;
    J[27] = -1.0*dwdx6;
    J[29] = -1.0*dwdx8;
    J[33] = -1.0*dwdx1;
    J[34] = -1.0*dwdx4;
    J[36] = -1.0*dwdx7;
    J[40] = -1.0*dwdx0;
    J[41] = 1.0*dwdx2;
    J[43] = 1.0*dwdx6;
    J[45] = 1.0*dwdx8 - 1.0*dwdx9;
    J[46] = -1.0*dwdx10;
    J[48] = 1.0*dwdx0;
    J[49] = -1.0*dwdx3;
    J[53] = 1.0*dwdx9;
    J[54] = 1.0*dwdx10 - 1.0*dwdx11;
    J[55] = -1.0*dwdx12;
    J[57] = 1.0*dwdx3;
    J[62] = 1.0*dwdx11;
    J[63] = 1.0*dwdx12;
}