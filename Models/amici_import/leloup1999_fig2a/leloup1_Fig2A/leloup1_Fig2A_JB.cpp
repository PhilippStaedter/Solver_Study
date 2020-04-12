#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_leloup1_Fig2A(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0 - 1.0*dwdx1 + 1.0*dwdx2;
    JB[1] = -1.0*dwdx2;
    JB[6] = 1.0*dwdx1;
    JB[9] = 1.0*dwdx1;
    JB[10] = 1.0*dwdx5;
    JB[11] = -1.0*dwdx5 + 1.0*dwdx6;
    JB[12] = -1.0*dwdx3;
    JB[13] = -1.0*dwdx4;
    JB[22] = 1.0*dwdx7;
    JB[24] = -1.0*dwdx8;
    JB[33] = 1.0*dwdx9;
    JB[37] = -1.0*dwdx10;
    JB[44] = 1.0*dwdx11 + 1.0*dwdx12;
    JB[45] = -1.0*dwdx12;
    JB[54] = -1.0*dwdx14;
    JB[55] = 1.0*dwdx13 + 1.0*dwdx14 + 1.0*dwdx15;
    JB[56] = -1.0*dwdx15;
    JB[60] = -1.0*dwdx18;
    JB[65] = -1.0*dwdx17;
    JB[66] = 1.0*dwdx16 + 1.0*dwdx17 + 1.0*dwdx18;
    JB[69] = 1.0*dwdx18;
    JB[77] = 1.0*dwdx19 + 1.0*dwdx20;
    JB[78] = -1.0*dwdx20;
    JB[87] = -1.0*dwdx22;
    JB[88] = 1.0*dwdx21 + 1.0*dwdx22 + 1.0*dwdx23;
    JB[89] = -1.0*dwdx23;
    JB[90] = -1.0*dwdx24;
    JB[96] = 1.0*dwdx24;
    JB[98] = -1.0*dwdx26;
    JB[99] = 1.0*dwdx24 + 1.0*dwdx25 + 1.0*dwdx26;
}