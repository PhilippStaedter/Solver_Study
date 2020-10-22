#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_leloup2_Fig2A(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx1 + 1.0*dwdx2;
    JB[1] = -1.0*dwdx0;
    JB[17] = 1.0*dwdx3 + 1.0*dwdx4 + 1.0*dwdx5;
    JB[24] = -1.0*dwdx3;
    JB[25] = -1.0*dwdx5;
    JB[34] = -1.0*dwdx6 + 1.0*dwdx7 + 1.0*dwdx8 + 1.0*dwdx9;
    JB[35] = -1.0*dwdx7;
    JB[36] = -1.0*dwdx8;
    JB[42] = 1.0*dwdx6;
    JB[46] = 1.0*dwdx6;
    JB[50] = -1.0*dwdx12;
    JB[51] = 1.0*dwdx10 + 1.0*dwdx11 + 1.0*dwdx12;
    JB[66] = 1.0*dwdx13;
    JB[68] = -1.0*dwdx13 + 1.0*dwdx14 + 1.0*dwdx15 + 1.0*dwdx16;
    JB[70] = -1.0*dwdx14;
    JB[71] = -1.0*dwdx15;
    JB[73] = 1.0*dwdx15;
    JB[85] = 1.0*dwdx17 + 1.0*dwdx18 + 1.0*dwdx19;
    JB[89] = -1.0*dwdx19;
    JB[100] = -1.0*dwdx22;
    JB[102] = 1.0*dwdx20 + 1.0*dwdx21 + 1.0*dwdx22;
    JB[116] = 1.0*dwdx23;
    JB[119] = -1.0*dwdx23 + 1.0*dwdx24 + 1.0*dwdx25;
    JB[121] = 1.0*dwdx23;
    JB[129] = -1.0*dwdx28;
    JB[136] = 1.0*dwdx26 + 1.0*dwdx27 + 1.0*dwdx28;
    JB[144] = -1.0*dwdx29;
    JB[145] = 1.0*dwdx33;
    JB[148] = 1.0*dwdx31;
    JB[149] = -1.0*dwdx32;
    JB[151] = -1.0*dwdx31;
    JB[153] = 1.0*dwdx31 + 1.0*dwdx32 - 1.0*dwdx33 + 1.0*dwdx34;
    JB[155] = -1.0*dwdx30;
    JB[157] = -1.0*dwdx35;
    JB[162] = -1.0*dwdx36;
    JB[170] = 1.0*dwdx36 + 1.0*dwdx37 + 1.0*dwdx38;
    JB[172] = -1.0*dwdx38;
    JB[174] = 1.0*dwdx36;
    JB[186] = -1.0*dwdx40;
    JB[187] = 1.0*dwdx39 + 1.0*dwdx41;
    JB[202] = -1.0*dwdx43;
    JB[204] = 1.0*dwdx42 + 1.0*dwdx43 + 1.0*dwdx44;
    JB[221] = 1.0*dwdx46 + 1.0*dwdx47;
    JB[222] = -1.0*dwdx45;
    JB[226] = -1.0*dwdx49;
    JB[234] = 1.0*dwdx49;
    JB[238] = 1.0*dwdx48 + 1.0*dwdx49 + 1.0*dwdx50;
    JB[239] = -1.0*dwdx48;
    JB[254] = -1.0*dwdx53;
    JB[255] = 1.0*dwdx51 + 1.0*dwdx52 + 1.0*dwdx53;
}