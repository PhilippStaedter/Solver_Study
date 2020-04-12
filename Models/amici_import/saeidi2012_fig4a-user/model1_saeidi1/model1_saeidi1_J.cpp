#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model1_saeidi1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0;
    J[4] = -1.0*dwdx6;
    J[11] = -1.0*dwdx1;
    J[12] = -1.0*dwdx3;
    J[17] = -1.0*dwdx11;
    J[21] = 1.0*dwdx1;
    J[22] = -1.0*dwdx2 + 1.0*dwdx3;
    J[27] = 1.0*dwdx11;
    J[28] = -1.0*dwdx12;
    J[33] = -1.0*dwdx4 - 1.0*dwdx5;
    J[34] = 1.0*dwdx7;
    J[36] = -1.0*dwdx8;
    J[37] = -1.0*dwdx10;
    J[40] = 1.0*dwdx0;
    J[44] = 1.0*dwdx6 - 1.0*dwdx7;
    J[53] = 1.0*dwdx4;
    J[63] = -1.0*dwdx5;
    J[66] = -1.0*dwdx8 - 1.0*dwdx9;
    J[67] = -1.0*dwdx10;
    J[71] = -1.0*dwdx1;
    J[72] = -1.0*dwdx3;
    J[73] = 1.0*dwdx5;
    J[76] = 1.0*dwdx8;
    J[77] = 1.0*dwdx10 - 1.0*dwdx11;
    J[82] = 1.0*dwdx2;
    J[88] = 1.0*dwdx12;
    J[96] = 1.0*dwdx9;
}