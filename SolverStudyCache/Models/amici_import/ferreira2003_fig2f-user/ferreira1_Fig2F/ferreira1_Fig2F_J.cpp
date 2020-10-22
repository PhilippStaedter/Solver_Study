#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_ferreira1_Fig2F(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 - 1.0*dwdx1;
    J[5] = 1.0*dwdx9;
    J[6] = 1.0*dwdx1;
    J[9] = 1.0*dwdx4;
    J[10] = 1.0*dwdx7;
    J[11] = 1.0*dwdx10;
    J[14] = -1.0*dwdx2 - 1.0*dwdx3;
    J[16] = -1.0*dwdx6;
    J[17] = 1.0*dwdx8;
    J[20] = 1.0*dwdx3;
    J[21] = -1.0*dwdx4 - 1.0*dwdx5;
    J[22] = -1.0*dwdx7;
    J[23] = 1.0*dwdx12;
    J[26] = -1.0*dwdx2;
    J[27] = -1.0*dwdx4;
    J[28] = -1.0*dwdx6 - 1.0*dwdx7;
    J[29] = 1.0*dwdx11 + 1.0*dwdx8;
    J[30] = 1.0*dwdx0;
    J[32] = 1.0*dwdx2;
    J[34] = 1.0*dwdx6;
    J[35] = -1.0*dwdx10 - 1.0*dwdx13 - 1.0*dwdx8 - 1.0*dwdx9;
}