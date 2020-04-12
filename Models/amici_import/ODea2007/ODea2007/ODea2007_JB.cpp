#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_ODea2007(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0;
    JB[1] = -1.0*dwdx1;
    JB[25] = 1.0*dwdx2 + 1.0*dwdx3 + 1.0*dwdx4 + 1.0*dwdx5;
    JB[26] = -1.0*dwdx5;
    JB[27] = -1.0*dwdx2;
    JB[28] = -1.0*dwdx3;
    JB[31] = 1.0*dwdx3;
    JB[32] = 1.0*dwdx2;
    JB[49] = 1.0*dwdx8;
    JB[50] = 1.0*dwdx6 + 1.0*dwdx7 - 1.0*dwdx8;
    JB[53] = -1.0*dwdx6;
    JB[57] = 1.0*dwdx6;
    JB[73] = 1.0*dwdx9;
    JB[75] = 1.0*dwdx10 + 1.0*dwdx11 - 1.0*dwdx9;
    JB[78] = -1.0*dwdx10;
    JB[79] = 1.0*dwdx10;
    JB[80] = -1.0*dwdx11 + 1.0*dwdx9;
    JB[97] = 1.0*dwdx12;
    JB[100] = -1.0*dwdx12 + 1.0*dwdx13 + 1.0*dwdx14;
    JB[102] = -1.0*dwdx13;
    JB[103] = 1.0*dwdx12 - 1.0*dwdx14;
    JB[104] = 1.0*dwdx13;
    JB[122] = 1.0*dwdx15;
    JB[124] = -1.0*dwdx17;
    JB[125] = -1.0*dwdx15 + 1.0*dwdx16 + 1.0*dwdx17;
    JB[129] = 1.0*dwdx15 - 1.0*dwdx16;
    JB[147] = 1.0*dwdx19;
    JB[148] = 1.0*dwdx18;
    JB[150] = -1.0*dwdx18 - 1.0*dwdx19 + 1.0*dwdx20;
    JB[151] = 1.0*dwdx19 - 1.0*dwdx20;
    JB[152] = 1.0*dwdx18 - 1.0*dwdx20;
    JB[169] = 1.0*dwdx21;
    JB[171] = 1.0*dwdx22;
    JB[172] = -1.0*dwdx21;
    JB[174] = -1.0*dwdx22;
    JB[175] = 1.0*dwdx21 + 1.0*dwdx22 + 1.0*dwdx23 + 1.0*dwdx24 + 1.0*dwdx25 + 1.0*dwdx26 + 1.0*dwdx27;
    JB[177] = -1.0*dwdx25;
    JB[178] = 1.0*dwdx23;
    JB[179] = -1.0*dwdx23;
    JB[181] = -1.0*dwdx24;
    JB[183] = 1.0*dwdx24;
    JB[186] = 1.0*dwdx26;
    JB[188] = -1.0*dwdx26;
    JB[190] = -1.0*dwdx27;
    JB[191] = 1.0*dwdx27;
    JB[193] = 1.0*dwdx28;
    JB[195] = -1.0*dwdx28;
    JB[196] = 1.0*dwdx29;
    JB[198] = -1.0*dwdx29;
    JB[200] = 1.0*dwdx28 + 1.0*dwdx29 + 1.0*dwdx30 + 1.0*dwdx31 + 1.0*dwdx32 + 1.0*dwdx33 + 1.0*dwdx34;
    JB[202] = -1.0*dwdx31;
    JB[203] = -1.0*dwdx30;
    JB[205] = 1.0*dwdx30;
    JB[207] = 1.0*dwdx31;
    JB[210] = 1.0*dwdx32;
    JB[212] = 1.0*dwdx33;
    JB[214] = -1.0*dwdx33;
    JB[215] = -1.0*dwdx32;
    JB[216] = -1.0*dwdx39;
    JB[218] = 1.0*dwdx35;
    JB[221] = -1.0*dwdx35;
    JB[223] = 1.0*dwdx37;
    JB[225] = 1.0*dwdx35 + 1.0*dwdx36 - 1.0*dwdx37 + 1.0*dwdx38;
    JB[228] = -1.0*dwdx36;
    JB[230] = 1.0*dwdx36;
    JB[235] = 1.0*dwdx38;
    JB[237] = -1.0*dwdx38;
    JB[247] = 1.0*dwdx41;
    JB[248] = -1.0*dwdx40 + 1.0*dwdx42;
    JB[250] = 1.0*dwdx40 + 1.0*dwdx41 - 1.0*dwdx42;
    JB[251] = -1.0*dwdx41;
    JB[255] = 1.0*dwdx42;
    JB[271] = -1.0*dwdx43 + 1.0*dwdx44;
    JB[272] = -1.0*dwdx43 + 1.0*dwdx45;
    JB[274] = 1.0*dwdx44;
    JB[275] = 1.0*dwdx43 - 1.0*dwdx44 - 1.0*dwdx45;
    JB[277] = 1.0*dwdx45;
    JB[297] = -1.0*dwdx47 + 1.0*dwdx48;
    JB[300] = 1.0*dwdx46 + 1.0*dwdx47 - 1.0*dwdx48;
    JB[301] = -1.0*dwdx46;
    JB[302] = 1.0*dwdx48;
    JB[319] = -1.0*dwdx49 + 1.0*dwdx51;
    JB[320] = 1.0*dwdx50;
    JB[323] = -1.0*dwdx50;
    JB[325] = 1.0*dwdx49 + 1.0*dwdx50 - 1.0*dwdx51;
    JB[327] = 1.0*dwdx51;
    JB[345] = 1.0*dwdx54;
    JB[348] = -1.0*dwdx54;
    JB[350] = -1.0*dwdx52 + 1.0*dwdx53 + 1.0*dwdx54;
    JB[351] = 1.0*dwdx52;
    JB[367] = 1.0*dwdx57;
    JB[368] = 1.0*dwdx58;
    JB[370] = -1.0*dwdx58;
    JB[373] = -1.0*dwdx57;
    JB[374] = -1.0*dwdx55;
    JB[375] = 1.0*dwdx55 + 1.0*dwdx56 + 1.0*dwdx57 + 1.0*dwdx58;
    JB[399] = -1.0*dwdx59;
    JB[400] = 1.0*dwdx60;
    JB[425] = 1.0*dwdx61;
    JB[426] = -1.0*dwdx62;
    JB[439] = 1.0*dwdx64;
    JB[440] = 1.0*dwdx63;
    JB[450] = 1.0*dwdx63 + 1.0*dwdx64 + 1.0*dwdx65 + 1.0*dwdx66;
    JB[451] = -1.0*dwdx66;
    JB[452] = -1.0*dwdx64;
    JB[455] = -1.0*dwdx63;
    JB[465] = 1.0*dwdx67;
    JB[474] = 1.0*dwdx69;
    JB[475] = 1.0*dwdx67 + 1.0*dwdx68 - 1.0*dwdx69;
    JB[477] = -1.0*dwdx67;
    JB[487] = 1.0*dwdx70 - 1.0*dwdx72;
    JB[488] = 1.0*dwdx71;
    JB[498] = 1.0*dwdx70;
    JB[500] = -1.0*dwdx70 + 1.0*dwdx71 + 1.0*dwdx72;
    JB[502] = -1.0*dwdx71;
    JB[513] = 1.0*dwdx73 - 1.0*dwdx74;
    JB[523] = 1.0*dwdx73;
    JB[524] = -1.0*dwdx75;
    JB[525] = -1.0*dwdx73 + 1.0*dwdx74 + 1.0*dwdx75;
    JB[535] = 1.0*dwdx77 - 1.0*dwdx78;
    JB[536] = 1.0*dwdx76 - 1.0*dwdx78;
    JB[548] = 1.0*dwdx76;
    JB[550] = -1.0*dwdx76 - 1.0*dwdx77 + 1.0*dwdx78;
    JB[551] = 1.0*dwdx77;
    JB[559] = 1.0*dwdx80;
    JB[560] = 1.0*dwdx79 - 1.0*dwdx81;
    JB[570] = 1.0*dwdx79;
    JB[574] = -1.0*dwdx80;
    JB[575] = -1.0*dwdx79 + 1.0*dwdx80 + 1.0*dwdx81;
}