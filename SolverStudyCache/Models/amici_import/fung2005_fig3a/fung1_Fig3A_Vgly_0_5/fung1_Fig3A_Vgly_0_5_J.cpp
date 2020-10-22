#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_fung1_Fig3A_Vgly_0_5(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 - 1.0*dwdx1;
    J[2] = 1.0*dwdx6;
    J[6] = 1.0*dwdx14;
    J[7] = -1.0*dwdx16;
    J[8] = 1.0*dwdx0;
    J[9] = -1.0*dwdx4;
    J[14] = -1.0*dwdx13;
    J[15] = 1.0*dwdx16;
    J[17] = 1.0*dwdx2;
    J[18] = -1.0*dwdx5;
    J[27] = 1.0*dwdx7 - 1.0*dwdx8;
    J[28] = -1.0*dwdx9;
    J[30] = 1.0*dwdx12;
    J[41] = 1.0*dwdx3;
    J[45] = -1.0*dwdx11;
    J[49] = 1.0*dwdx4;
    J[50] = -1.0*dwdx6;
    J[51] = -1.0*dwdx7;
    J[54] = -1.0*dwdx12 + 1.0*dwdx13 - 1.0*dwdx14;
    J[61] = 1.0*dwdx10;
    J[63] = -1.0*dwdx15;
}