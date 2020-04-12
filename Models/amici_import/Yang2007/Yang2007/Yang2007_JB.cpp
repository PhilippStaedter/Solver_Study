#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_Yang2007(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = -1.0*dwdx0 + 1.0*dwdx1 + 1.0*dwdx2 + 1.0*dwdx3 + 1.0*dwdx5 + 1.0*dwdx6;
    JB[1] = -1.0*dwdx5;
    JB[11] = -1.0*dwdx1;
    JB[19] = -1.0*dwdx2;
    JB[21] = -1.0*dwdx3 + 1.0*dwdx4;
    JB[22] = -1.0*dwdx4;
    JB[25] = 1.0*dwdx7;
    JB[26] = -1.0*dwdx7 + 1.0*dwdx8 + 1.0*dwdx9;
    JB[27] = -1.0*dwdx8;
    JB[28] = -1.0*dwdx9;
    JB[38] = 1.0*dwdx10;
    JB[50] = -1.0*dwdx11 + 1.0*dwdx12;
    JB[51] = -1.0*dwdx12 + 1.0*dwdx13 + 1.0*dwdx14;
    JB[52] = -1.0*dwdx13 + 1.0*dwdx16;
    JB[53] = -1.0*dwdx14;
    JB[54] = 1.0*dwdx15;
    JB[55] = -1.0*dwdx15;
    JB[76] = 1.0*dwdx17;
    JB[78] = -1.0*dwdx17 + 1.0*dwdx18 + 1.0*dwdx20;
    JB[79] = -1.0*dwdx18;
    JB[88] = 1.0*dwdx21;
    JB[89] = 1.0*dwdx19;
    JB[100] = -1.0*dwdx22;
    JB[103] = 1.0*dwdx23;
    JB[104] = -1.0*dwdx23 + 1.0*dwdx24 + 1.0*dwdx26;
    JB[105] = -1.0*dwdx24;
    JB[113] = -1.0*dwdx25;
    JB[129] = 1.0*dwdx27;
    JB[130] = -1.0*dwdx27;
    JB[150] = -1.0*dwdx28;
    JB[175] = 1.0*dwdx29;
    JB[182] = 1.0*dwdx30;
    JB[186] = -1.0*dwdx29;
    JB[200] = 1.0*dwdx31;
    JB[208] = 1.0*dwdx32;
    JB[219] = -1.0*dwdx31;
    JB[225] = 1.0*dwdx33;
    JB[246] = -1.0*dwdx33;
    JB[271] = 1.0*dwdx34;
    JB[272] = -1.0*dwdx34;
    JB[275] = -1.0*dwdx35 + 1.0*dwdx36;
    JB[283] = 1.0*dwdx38;
    JB[286] = -1.0*dwdx36 + 1.0*dwdx37 + 1.0*dwdx40;
    JB[287] = 1.0*dwdx39;
    JB[288] = 1.0*dwdx41;
    JB[293] = -1.0*dwdx37;
    JB[312] = 1.0*dwdx43 + 1.0*dwdx44;
    JB[321] = 1.0*dwdx42;
    JB[323] = -1.0*dwdx42;
    JB[325] = 1.0*dwdx45;
    JB[326] = -1.0*dwdx45 + 1.0*dwdx46;
    JB[328] = -1.0*dwdx46;
    JB[338] = -1.0*dwdx47 + 1.0*dwdx48 + 1.0*dwdx49 + 1.0*dwdx50;
    JB[353] = 1.0*dwdx51;
    JB[354] = -1.0*dwdx51;
    JB[364] = 1.0*dwdx52;
    JB[379] = 1.0*dwdx53;
    JB[380] = -1.0*dwdx53;
    JB[401] = 1.0*dwdx56;
    JB[402] = -1.0*dwdx56;
    JB[411] = 1.0*dwdx54;
    JB[418] = -1.0*dwdx54;
    JB[419] = 1.0*dwdx55;
    JB[420] = -1.0*dwdx55;
    JB[450] = 1.0*dwdx58 + 1.0*dwdx60;
    JB[451] = -1.0*dwdx60 + 1.0*dwdx61;
    JB[453] = -1.0*dwdx61;
    JB[461] = 1.0*dwdx57;
    JB[468] = -1.0*dwdx57 + 1.0*dwdx62;
    JB[469] = -1.0*dwdx58;
    JB[471] = 1.0*dwdx59;
    JB[472] = -1.0*dwdx59;
    JB[475] = -1.0*dwdx63 + 1.0*dwdx64;
    JB[494] = -1.0*dwdx64 + 1.0*dwdx65;
    JB[495] = -1.0*dwdx65;
    JB[500] = 1.0*dwdx67;
    JB[501] = -1.0*dwdx67 + 1.0*dwdx68;
    JB[503] = -1.0*dwdx68;
    JB[504] = 1.0*dwdx69;
    JB[505] = -1.0*dwdx69;
    JB[519] = 1.0*dwdx66;
    JB[520] = -1.0*dwdx66;
    JB[525] = 1.0*dwdx70;
    JB[537] = 1.0*dwdx73;
    JB[546] = -1.0*dwdx70 + 1.0*dwdx71 + 1.0*dwdx72;
    JB[547] = -1.0*dwdx71;
    JB[548] = -1.0*dwdx72;
    JB[550] = 1.0*dwdx74 + 1.0*dwdx76;
    JB[551] = -1.0*dwdx76 + 1.0*dwdx77;
    JB[553] = -1.0*dwdx77;
    JB[557] = -1.0*dwdx78;
    JB[571] = -1.0*dwdx74 + 1.0*dwdx75;
    JB[572] = -1.0*dwdx75;
    JB[596] = 1.0*dwdx79;
    JB[598] = -1.0*dwdx79 + 1.0*dwdx81;
    JB[599] = -1.0*dwdx80;
    JB[624] = 1.0*dwdx82;
}