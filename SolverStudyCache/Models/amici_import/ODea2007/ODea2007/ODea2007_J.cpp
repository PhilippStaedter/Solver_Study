#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_ODea2007(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0;
    J[9] = 1.0*dwdx39;
    J[24] = 1.0*dwdx1;
    J[25] = -1.0*dwdx2 - 1.0*dwdx3 - 1.0*dwdx4 - 1.0*dwdx5;
    J[26] = -1.0*dwdx8;
    J[27] = -1.0*dwdx9;
    J[28] = -1.0*dwdx12;
    J[31] = -1.0*dwdx21;
    J[32] = -1.0*dwdx28;
    J[49] = 1.0*dwdx5;
    J[50] = -1.0*dwdx6 - 1.0*dwdx7 + 1.0*dwdx8;
    J[53] = -1.0*dwdx15;
    J[57] = -1.0*dwdx35;
    J[73] = 1.0*dwdx2;
    J[75] = -1.0*dwdx10 - 1.0*dwdx11 + 1.0*dwdx9;
    J[78] = -1.0*dwdx19;
    J[79] = -1.0*dwdx22;
    J[80] = 1.0*dwdx28;
    J[97] = 1.0*dwdx3;
    J[100] = 1.0*dwdx12 - 1.0*dwdx13 - 1.0*dwdx14;
    J[101] = 1.0*dwdx17;
    J[102] = -1.0*dwdx18;
    J[103] = 1.0*dwdx21;
    J[104] = -1.0*dwdx29;
    J[122] = 1.0*dwdx6;
    J[125] = 1.0*dwdx15 - 1.0*dwdx16 - 1.0*dwdx17;
    J[129] = 1.0*dwdx35;
    J[147] = 1.0*dwdx10;
    J[148] = 1.0*dwdx13;
    J[150] = 1.0*dwdx18 + 1.0*dwdx19 - 1.0*dwdx20;
    J[151] = 1.0*dwdx22;
    J[152] = 1.0*dwdx29;
    J[169] = -1.0*dwdx3;
    J[171] = -1.0*dwdx10;
    J[172] = -1.0*dwdx12 + 1.0*dwdx14;
    J[174] = -1.0*dwdx19 + 1.0*dwdx20;
    J[175] = -1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24 - 1.0*dwdx25 - 1.0*dwdx26 - 1.0*dwdx27;
    J[177] = -1.0*dwdx37;
    J[178] = -1.0*dwdx41;
    J[179] = 1.0*dwdx43 - 1.0*dwdx44;
    J[181] = 1.0*dwdx49 - 1.0*dwdx51;
    J[183] = -1.0*dwdx57;
    J[186] = -1.0*dwdx64;
    J[188] = -1.0*dwdx70 + 1.0*dwdx72;
    J[190] = -1.0*dwdx77 + 1.0*dwdx78;
    J[191] = -1.0*dwdx80;
    J[193] = -1.0*dwdx2;
    J[195] = 1.0*dwdx11 - 1.0*dwdx9;
    J[196] = -1.0*dwdx13;
    J[198] = -1.0*dwdx18 + 1.0*dwdx20;
    J[200] = -1.0*dwdx28 - 1.0*dwdx29 - 1.0*dwdx30 - 1.0*dwdx31 - 1.0*dwdx32 - 1.0*dwdx33 - 1.0*dwdx34;
    J[202] = 1.0*dwdx40 - 1.0*dwdx42;
    J[203] = 1.0*dwdx43 - 1.0*dwdx45;
    J[205] = -1.0*dwdx50;
    J[207] = -1.0*dwdx58;
    J[210] = -1.0*dwdx63;
    J[212] = -1.0*dwdx71;
    J[214] = -1.0*dwdx76 + 1.0*dwdx78;
    J[215] = -1.0*dwdx79 + 1.0*dwdx81;
    J[218] = -1.0*dwdx6;
    J[221] = -1.0*dwdx15 + 1.0*dwdx16;
    J[223] = 1.0*dwdx25;
    J[225] = -1.0*dwdx35 - 1.0*dwdx36 + 1.0*dwdx37 - 1.0*dwdx38;
    J[228] = 1.0*dwdx47 - 1.0*dwdx48;
    J[230] = -1.0*dwdx54;
    J[235] = -1.0*dwdx67;
    J[237] = -1.0*dwdx73 + 1.0*dwdx74;
    J[247] = -1.0*dwdx23;
    J[248] = 1.0*dwdx31;
    J[250] = -1.0*dwdx40 - 1.0*dwdx41 + 1.0*dwdx42;
    J[251] = -1.0*dwdx44;
    J[255] = 1.0*dwdx58;
    J[271] = 1.0*dwdx23;
    J[272] = 1.0*dwdx30;
    J[274] = 1.0*dwdx41;
    J[275] = -1.0*dwdx43 + 1.0*dwdx44 + 1.0*dwdx45;
    J[277] = 1.0*dwdx50;
    J[297] = 1.0*dwdx36;
    J[300] = -1.0*dwdx46 - 1.0*dwdx47 + 1.0*dwdx48;
    J[302] = 1.0*dwdx54;
    J[319] = 1.0*dwdx24;
    J[320] = -1.0*dwdx30;
    J[323] = -1.0*dwdx45;
    J[324] = 1.0*dwdx46;
    J[325] = -1.0*dwdx49 - 1.0*dwdx50 + 1.0*dwdx51;
    J[327] = 1.0*dwdx57;
    J[345] = -1.0*dwdx36;
    J[348] = -1.0*dwdx48;
    J[350] = 1.0*dwdx52 - 1.0*dwdx53 - 1.0*dwdx54;
    J[351] = 1.0*dwdx55;
    J[367] = -1.0*dwdx24;
    J[368] = -1.0*dwdx31;
    J[370] = -1.0*dwdx42;
    J[373] = -1.0*dwdx51;
    J[374] = -1.0*dwdx52;
    J[375] = -1.0*dwdx55 - 1.0*dwdx56 - 1.0*dwdx57 - 1.0*dwdx58;
    J[376] = 1.0*dwdx59;
    J[400] = -1.0*dwdx60;
    J[425] = -1.0*dwdx61;
    J[439] = -1.0*dwdx26;
    J[440] = -1.0*dwdx32;
    J[449] = 1.0*dwdx62;
    J[450] = -1.0*dwdx63 - 1.0*dwdx64 - 1.0*dwdx65 - 1.0*dwdx66;
    J[451] = -1.0*dwdx69;
    J[452] = -1.0*dwdx70;
    J[455] = -1.0*dwdx79;
    J[465] = -1.0*dwdx38;
    J[474] = 1.0*dwdx66;
    J[475] = -1.0*dwdx67 - 1.0*dwdx68 + 1.0*dwdx69;
    J[477] = -1.0*dwdx73;
    J[487] = 1.0*dwdx26;
    J[488] = -1.0*dwdx33;
    J[498] = 1.0*dwdx64;
    J[500] = 1.0*dwdx70 - 1.0*dwdx71 - 1.0*dwdx72;
    J[501] = 1.0*dwdx75;
    J[502] = -1.0*dwdx76;
    J[513] = 1.0*dwdx38;
    J[523] = 1.0*dwdx67;
    J[525] = 1.0*dwdx73 - 1.0*dwdx74 - 1.0*dwdx75;
    J[535] = 1.0*dwdx27;
    J[536] = 1.0*dwdx33;
    J[548] = 1.0*dwdx71;
    J[550] = 1.0*dwdx76 + 1.0*dwdx77 - 1.0*dwdx78;
    J[551] = 1.0*dwdx80;
    J[559] = -1.0*dwdx27;
    J[560] = 1.0*dwdx32;
    J[570] = 1.0*dwdx63;
    J[574] = -1.0*dwdx77;
    J[575] = 1.0*dwdx79 - 1.0*dwdx80 - 1.0*dwdx81;
}