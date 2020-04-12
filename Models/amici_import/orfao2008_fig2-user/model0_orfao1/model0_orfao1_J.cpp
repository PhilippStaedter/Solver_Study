#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_orfao1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 - 1.0*dwdx1;
    J[5] = -1.0*dwdx7;
    J[10] = -1.0*dwdx14;
    J[12] = 1.0*dwdx0 + 1.0*dwdx1;
    J[13] = -1.0*dwdx3 - 1.0*dwdx4;
    J[17] = 1.0*dwdx7;
    J[22] = 1.0*dwdx14;
    J[25] = 1.0*dwdx4;
    J[37] = 1.0*dwdx3;
    J[52] = -1.0*dwdx5;
    J[53] = -1.0*dwdx6;
    J[56] = -1.0*dwdx10;
    J[58] = -1.0*dwdx13;
    J[64] = 1.0*dwdx5;
    J[65] = 1.0*dwdx6;
    J[68] = 1.0*dwdx10;
    J[70] = 1.0*dwdx13;
    J[85] = -1.0*dwdx2;
    J[91] = -1.0*dwdx9;
    J[97] = 1.0*dwdx2;
    J[100] = -1.0*dwdx5;
    J[101] = -1.0*dwdx6;
    J[103] = 1.0*dwdx9;
    J[104] = -1.0*dwdx10;
    J[106] = -1.0*dwdx13;
    J[114] = -1.0*dwdx8;
    J[117] = -1.0*dwdx11;
    J[124] = -1.0*dwdx5;
    J[125] = -1.0*dwdx6;
    J[126] = 1.0*dwdx8;
    J[128] = -1.0*dwdx10;
    J[129] = 1.0*dwdx11;
    J[130] = -1.0*dwdx12 - 1.0*dwdx13;
    J[142] = 1.0*dwdx12;
}