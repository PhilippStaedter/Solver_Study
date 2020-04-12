#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_nyman3(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0;
    J[3] = 1.0*dwdx4;
    J[5] = 1.0*dwdx7;
    J[6] = 1.0*dwdx10;
    J[10] = -1.0*dwdx1;
    J[11] = 1.0*dwdx3;
    J[13] = -1.0*dwdx6;
    J[15] = -1.0*dwdx11;
    J[17] = -1.0*dwdx14;
    J[19] = 1.0*dwdx1;
    J[20] = -1.0*dwdx3;
    J[22] = 1.0*dwdx6;
    J[24] = 1.0*dwdx11;
    J[26] = 1.0*dwdx14;
    J[30] = -1.0*dwdx4;
    J[31] = 1.0*dwdx5;
    J[35] = 1.0*dwdx13;
    J[40] = -1.0*dwdx5;
    J[42] = 1.0*dwdx9;
    J[44] = -1.0*dwdx13;
    J[45] = 1.0*dwdx0;
    J[50] = -1.0*dwdx7 - 1.0*dwdx8;
    J[59] = 1.0*dwdx8;
    J[60] = -1.0*dwdx10 - 1.0*dwdx9;
    J[65] = -1.0*dwdx2;
    J[70] = -1.0*dwdx12;
    J[71] = 1.0*dwdx15;
    J[74] = 1.0*dwdx2;
    J[79] = 1.0*dwdx12;
    J[80] = -1.0*dwdx15;
}