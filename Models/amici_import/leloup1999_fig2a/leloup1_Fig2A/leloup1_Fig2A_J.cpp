#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_leloup1_Fig2A(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 + 1.0*dwdx1 - 1.0*dwdx2;
    J[1] = -1.0*dwdx5;
    J[6] = 1.0*dwdx18;
    J[9] = 1.0*dwdx24;
    J[10] = 1.0*dwdx2;
    J[11] = 1.0*dwdx5 - 1.0*dwdx6;
    J[21] = 1.0*dwdx3;
    J[22] = -1.0*dwdx7;
    J[31] = 1.0*dwdx4;
    J[33] = -1.0*dwdx9;
    J[42] = 1.0*dwdx8;
    J[44] = -1.0*dwdx11 - 1.0*dwdx12;
    J[45] = 1.0*dwdx14;
    J[54] = 1.0*dwdx12;
    J[55] = -1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx15;
    J[56] = 1.0*dwdx17;
    J[60] = -1.0*dwdx1;
    J[65] = 1.0*dwdx15;
    J[66] = -1.0*dwdx16 - 1.0*dwdx17 - 1.0*dwdx18;
    J[69] = -1.0*dwdx24;
    J[73] = 1.0*dwdx10;
    J[77] = -1.0*dwdx19 - 1.0*dwdx20;
    J[78] = 1.0*dwdx22;
    J[87] = 1.0*dwdx20;
    J[88] = -1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23;
    J[89] = 1.0*dwdx26;
    J[90] = -1.0*dwdx1;
    J[96] = -1.0*dwdx18;
    J[98] = 1.0*dwdx23;
    J[99] = -1.0*dwdx24 - 1.0*dwdx25 - 1.0*dwdx26;
}