#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_xie1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2;
    J[1] = 1.0*dwdx8;
    J[2] = 1.0*dwdx9;
    J[3] = 1.0*dwdx11;
    J[6] = -1.0*dwdx16;
    J[25] = 1.0*dwdx0;
    J[26] = -1.0*dwdx7 - 1.0*dwdx8;
    J[31] = 1.0*dwdx16;
    J[50] = 1.0*dwdx2;
    J[52] = -1.0*dwdx10 - 1.0*dwdx9;
    J[53] = -1.0*dwdx11;
    J[59] = 1.0*dwdx23;
    J[104] = -1.0*dwdx12;
    J[111] = 1.0*dwdx27;
    J[130] = -1.0*dwdx14 - 1.0*dwdx15;
    J[131] = 1.0*dwdx18;
    J[132] = -1.0*dwdx19;
    J[138] = 1.0*dwdx30;
    J[150] = -1.0*dwdx0;
    J[151] = 1.0*dwdx8;
    J[155] = 1.0*dwdx14;
    J[156] = -1.0*dwdx16 - 1.0*dwdx17 - 1.0*dwdx18;
    J[157] = 1.0*dwdx19;
    J[180] = -1.0*dwdx14;
    J[181] = 1.0*dwdx18;
    J[182] = -1.0*dwdx19 - 1.0*dwdx20;
    J[196] = 1.0*dwdx53;
    J[208] = -1.0*dwdx21;
    J[223] = 1.0*dwdx56;
    J[234] = -1.0*dwdx24;
    J[235] = 1.0*dwdx25;
    J[244] = 1.0*dwdx44;
    J[245] = 1.0*dwdx48;
    J[286] = -1.0*dwdx26;
    J[287] = 1.0*dwdx28;
    J[290] = 1.0*dwdx32;
    J[338] = -1.0*dwdx29;
    J[339] = 1.0*dwdx31;
    J[341] = 1.0*dwdx35;
    J[375] = 1.0*dwdx5;
    J[390] = -1.0*dwdx33 + 1.0*dwdx34;
    J[400] = 1.0*dwdx3;
    J[416] = 1.0*dwdx36 - 1.0*dwdx37;
    J[425] = 1.0*dwdx6;
    J[442] = 1.0*dwdx39 - 1.0*dwdx40;
    J[450] = 1.0*dwdx4;
    J[468] = -1.0*dwdx42 + 1.0*dwdx43;
    J[479] = 1.0*dwdx13;
    J[494] = 1.0*dwdx46 - 1.0*dwdx47;
    J[495] = 1.0*dwdx51;
    J[508] = 1.0*dwdx22;
    J[519] = 1.0*dwdx45;
    J[520] = 1.0*dwdx49 - 1.0*dwdx50;
    J[542] = 1.0*dwdx38;
    J[546] = -1.0*dwdx52;
    J[547] = 1.0*dwdx54;
    J[593] = 1.0*dwdx41;
    J[598] = -1.0*dwdx55;
    J[599] = 1.0*dwdx57;
}