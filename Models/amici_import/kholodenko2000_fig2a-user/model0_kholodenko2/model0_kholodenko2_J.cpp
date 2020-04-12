#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_kholodenko2(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0;
    J[1] = 1.0*dwdx2;
    J[7] = -1.0*dwdx13;
    J[8] = 1.0*dwdx0;
    J[9] = -1.0*dwdx1 - 1.0*dwdx2;
    J[10] = 1.0*dwdx4;
    J[15] = 1.0*dwdx13 - 1.0*dwdx14;
    J[17] = 1.0*dwdx1;
    J[18] = -1.0*dwdx4;
    J[23] = 1.0*dwdx14;
    J[27] = -1.0*dwdx5;
    J[29] = -1.0*dwdx8;
    J[30] = 1.0*dwdx11;
    J[34] = -1.0*dwdx3;
    J[36] = -1.0*dwdx6;
    J[37] = 1.0*dwdx7;
    J[42] = 1.0*dwdx3;
    J[44] = 1.0*dwdx6;
    J[45] = -1.0*dwdx7;
    J[51] = 1.0*dwdx5;
    J[53] = 1.0*dwdx8 - 1.0*dwdx9;
    J[54] = -1.0*dwdx10 - 1.0*dwdx11;
    J[55] = 1.0*dwdx12;
    J[61] = 1.0*dwdx9;
    J[62] = 1.0*dwdx10;
    J[63] = -1.0*dwdx12;
}