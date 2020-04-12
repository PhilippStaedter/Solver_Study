#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_kowald1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3 - 1.0*dwdx4;
    J[1] = -1.0*dwdx6 - 1.0*dwdx7;
    J[2] = -1.0*dwdx11;
    J[3] = -1.0*dwdx16;
    J[7] = -1.0*dwdx21;
    J[9] = -1.0*dwdx1 + 1.0*dwdx2;
    J[10] = -1.0*dwdx6 + 1.0*dwdx7 + 1.0*dwdx8 - 1.0*dwdx9;
    J[16] = 1.0*dwdx21 + 1.0*dwdx22;
    J[18] = 1.0*dwdx2 - 1.0*dwdx4 + 1.0*dwdx5;
    J[19] = -1.0*dwdx10 + 1.0*dwdx7;
    J[20] = -1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx13;
    J[25] = 1.0*dwdx21;
    J[26] = -1.0*dwdx23;
    J[27] = -1.0*dwdx3;
    J[30] = -1.0*dwdx14 - 2.0*dwdx15 - 1.0*dwdx16;
    J[33] = 1.0*dwdx20;
    J[36] = 2.0*dwdx4;
    J[37] = 2.0*dwdx10;
    J[38] = 2.0*dwdx11 + 2.0*dwdx12;
    J[40] = -1.0*dwdx17 - 1.0*dwdx18;
    J[45] = 1.0*dwdx3;
    J[48] = 1.0*dwdx14 + 1.0*dwdx16;
    J[50] = -1.0*dwdx19;
    J[54] = 1.0*dwdx5;
    J[57] = 1.0*dwdx14;
    J[58] = 1.0*dwdx18;
    J[60] = -1.0*dwdx20;
}