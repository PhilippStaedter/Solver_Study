#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_leloup2_Fig2A(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx1 - 1.0*dwdx2;
    J[9] = 1.0*dwdx29;
    J[16] = 1.0*dwdx0;
    J[17] = -1.0*dwdx3 - 1.0*dwdx4 - 1.0*dwdx5;
    J[24] = 1.0*dwdx28;
    J[25] = -1.0*dwdx33;
    J[34] = 1.0*dwdx6 - 1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    J[35] = 1.0*dwdx12;
    J[36] = -1.0*dwdx13;
    J[42] = 1.0*dwdx36;
    J[46] = 1.0*dwdx49;
    J[50] = 1.0*dwdx7;
    J[51] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12;
    J[66] = 1.0*dwdx8;
    J[68] = 1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx15 - 1.0*dwdx16;
    J[70] = 1.0*dwdx22;
    J[71] = -1.0*dwdx23;
    J[73] = -1.0*dwdx31;
    J[85] = -1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19;
    J[89] = 1.0*dwdx32;
    J[100] = 1.0*dwdx14;
    J[102] = -1.0*dwdx20 - 1.0*dwdx21 - 1.0*dwdx22;
    J[116] = 1.0*dwdx15;
    J[119] = 1.0*dwdx23 - 1.0*dwdx24 - 1.0*dwdx25;
    J[121] = 1.0*dwdx31;
    J[129] = 1.0*dwdx3;
    J[136] = -1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28;
    J[145] = 1.0*dwdx5;
    J[148] = -1.0*dwdx15;
    J[149] = 1.0*dwdx19;
    J[151] = -1.0*dwdx23;
    J[153] = -1.0*dwdx31 - 1.0*dwdx32 + 1.0*dwdx33 - 1.0*dwdx34;
    J[162] = -1.0*dwdx6;
    J[170] = -1.0*dwdx36 - 1.0*dwdx37 - 1.0*dwdx38;
    J[171] = 1.0*dwdx40;
    J[172] = 1.0*dwdx43;
    J[174] = -1.0*dwdx49;
    J[185] = 1.0*dwdx30;
    J[187] = -1.0*dwdx39 - 1.0*dwdx41;
    J[202] = 1.0*dwdx38;
    J[204] = -1.0*dwdx42 - 1.0*dwdx43 - 1.0*dwdx44;
    J[217] = 1.0*dwdx35;
    J[221] = -1.0*dwdx46 - 1.0*dwdx47;
    J[226] = -1.0*dwdx6;
    J[234] = -1.0*dwdx36;
    J[237] = 1.0*dwdx45;
    J[238] = -1.0*dwdx48 - 1.0*dwdx49 - 1.0*dwdx50;
    J[239] = 1.0*dwdx53;
    J[254] = 1.0*dwdx48;
    J[255] = -1.0*dwdx51 - 1.0*dwdx52 - 1.0*dwdx53;
}