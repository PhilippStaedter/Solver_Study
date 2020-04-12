#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_xie1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0 + 1.0*dwdx1 + 1.0*dwdx2;
    JB[1] = -1.0*dwdx0;
    JB[2] = -1.0*dwdx2;
    JB[6] = 1.0*dwdx0;
    JB[15] = -1.0*dwdx5;
    JB[16] = -1.0*dwdx3;
    JB[17] = -1.0*dwdx6;
    JB[18] = -1.0*dwdx4;
    JB[25] = -1.0*dwdx8;
    JB[26] = 1.0*dwdx7 + 1.0*dwdx8;
    JB[31] = -1.0*dwdx8;
    JB[50] = -1.0*dwdx9;
    JB[52] = 1.0*dwdx10 + 1.0*dwdx9;
    JB[75] = -1.0*dwdx11;
    JB[77] = 1.0*dwdx11;
    JB[104] = 1.0*dwdx12;
    JB[119] = -1.0*dwdx13;
    JB[130] = 1.0*dwdx14 + 1.0*dwdx15;
    JB[131] = -1.0*dwdx14;
    JB[132] = 1.0*dwdx14;
    JB[150] = 1.0*dwdx16;
    JB[151] = -1.0*dwdx16;
    JB[155] = -1.0*dwdx18;
    JB[156] = 1.0*dwdx16 + 1.0*dwdx17 + 1.0*dwdx18;
    JB[157] = -1.0*dwdx18;
    JB[180] = 1.0*dwdx19;
    JB[181] = -1.0*dwdx19;
    JB[182] = 1.0*dwdx19 + 1.0*dwdx20;
    JB[208] = 1.0*dwdx21;
    JB[220] = -1.0*dwdx22;
    JB[227] = -1.0*dwdx23;
    JB[234] = 1.0*dwdx24;
    JB[259] = -1.0*dwdx25;
    JB[279] = -1.0*dwdx27;
    JB[286] = 1.0*dwdx26;
    JB[311] = -1.0*dwdx28;
    JB[330] = -1.0*dwdx30;
    JB[338] = 1.0*dwdx29;
    JB[363] = -1.0*dwdx31;
    JB[386] = -1.0*dwdx32;
    JB[390] = 1.0*dwdx33 - 1.0*dwdx34;
    JB[413] = -1.0*dwdx35;
    JB[416] = -1.0*dwdx36 + 1.0*dwdx37;
    JB[442] = -1.0*dwdx39 + 1.0*dwdx40;
    JB[446] = -1.0*dwdx38;
    JB[468] = 1.0*dwdx42 - 1.0*dwdx43;
    JB[473] = -1.0*dwdx41;
    JB[484] = -1.0*dwdx44;
    JB[494] = -1.0*dwdx46 + 1.0*dwdx47;
    JB[495] = -1.0*dwdx45;
    JB[509] = -1.0*dwdx48;
    JB[519] = -1.0*dwdx51;
    JB[520] = -1.0*dwdx49 + 1.0*dwdx50;
    JB[532] = -1.0*dwdx53;
    JB[546] = 1.0*dwdx52;
    JB[571] = -1.0*dwdx54;
    JB[583] = -1.0*dwdx56;
    JB[598] = 1.0*dwdx55;
    JB[623] = -1.0*dwdx57;
}