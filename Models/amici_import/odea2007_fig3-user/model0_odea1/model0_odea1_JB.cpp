#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_odea1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0 + 1.0*dwdx1 + 1.0*dwdx2 + 1.0*dwdx3 + 1.0*dwdx4 + 1.0*dwdx5 + 1.0*dwdx6;
    JB[1] = -1.0*dwdx4;
    JB[2] = -1.0*dwdx1;
    JB[3] = 1.0*dwdx1;
    JB[5] = 1.0*dwdx4;
    JB[8] = -1.0*dwdx5;
    JB[9] = -1.0*dwdx2;
    JB[10] = 1.0*dwdx2;
    JB[12] = 1.0*dwdx5;
    JB[15] = -1.0*dwdx6;
    JB[16] = -1.0*dwdx3;
    JB[17] = 1.0*dwdx3;
    JB[19] = 1.0*dwdx6;
    JB[24] = -1.0*dwdx7 + 1.0*dwdx9;
    JB[25] = 1.0*dwdx7 + 1.0*dwdx8 - 1.0*dwdx9;
    JB[26] = -1.0*dwdx8;
    JB[29] = 1.0*dwdx9;
    JB[46] = 1.0*dwdx8;
    JB[48] = -1.0*dwdx10 + 1.0*dwdx12;
    JB[49] = 1.0*dwdx11;
    JB[50] = 1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12;
    JB[51] = 1.0*dwdx12;
    JB[70] = -1.0*dwdx10 + 1.0*dwdx11;
    JB[72] = 1.0*dwdx14;
    JB[74] = -1.0*dwdx14;
    JB[75] = 1.0*dwdx13 + 1.0*dwdx14 - 1.0*dwdx15;
    JB[77] = 1.0*dwdx15;
    JB[94] = -1.0*dwdx13 + 1.0*dwdx15;
    JB[99] = -1.0*dwdx18;
    JB[100] = 1.0*dwdx16 - 1.0*dwdx17 + 1.0*dwdx18;
    JB[103] = 1.0*dwdx17;
    JB[119] = -1.0*dwdx16 + 1.0*dwdx17;
    JB[120] = 1.0*dwdx20;
    JB[121] = -1.0*dwdx20;
    JB[123] = -1.0*dwdx21;
    JB[125] = 1.0*dwdx19 + 1.0*dwdx20 + 1.0*dwdx21 + 1.0*dwdx22;
    JB[127] = -1.0*dwdx22;
    JB[142] = 1.0*dwdx21;
    JB[149] = -1.0*dwdx24;
    JB[150] = 1.0*dwdx23;
    JB[172] = -1.0*dwdx26;
    JB[173] = 1.0*dwdx27;
    JB[175] = 1.0*dwdx25 + 1.0*dwdx26 - 1.0*dwdx27;
    JB[191] = 1.0*dwdx26;
    JB[192] = -1.0*dwdx28 + 1.0*dwdx30;
    JB[200] = 1.0*dwdx28 + 1.0*dwdx29 - 1.0*dwdx30;
    JB[201] = -1.0*dwdx29;
    JB[204] = 1.0*dwdx30;
    JB[214] = 1.0*dwdx29;
    JB[216] = -1.0*dwdx31 + 1.0*dwdx33;
    JB[224] = 1.0*dwdx32;
    JB[225] = 1.0*dwdx31 - 1.0*dwdx32 - 1.0*dwdx33;
    JB[226] = 1.0*dwdx33;
    JB[238] = -1.0*dwdx31 + 1.0*dwdx32;
    JB[240] = 1.0*dwdx35;
    JB[249] = -1.0*dwdx35;
    JB[250] = 1.0*dwdx34 + 1.0*dwdx35 - 1.0*dwdx36;
    JB[252] = 1.0*dwdx36;
    JB[262] = -1.0*dwdx34 + 1.0*dwdx36;
    JB[274] = -1.0*dwdx39;
    JB[275] = 1.0*dwdx37 - 1.0*dwdx38 + 1.0*dwdx39;
    JB[278] = 1.0*dwdx38;
    JB[287] = -1.0*dwdx37 + 1.0*dwdx38;
    JB[288] = 1.0*dwdx41;
    JB[296] = -1.0*dwdx41;
    JB[298] = -1.0*dwdx42;
    JB[300] = 1.0*dwdx40 + 1.0*dwdx41 + 1.0*dwdx42 + 1.0*dwdx43;
    JB[302] = -1.0*dwdx43;
    JB[310] = 1.0*dwdx42;
    JB[324] = -1.0*dwdx45;
    JB[325] = 1.0*dwdx44;
    JB[347] = -1.0*dwdx47;
    JB[348] = 1.0*dwdx48;
    JB[350] = 1.0*dwdx46 + 1.0*dwdx47 - 1.0*dwdx48;
    JB[359] = 1.0*dwdx47;
    JB[360] = -1.0*dwdx49 + 1.0*dwdx51;
    JB[375] = 1.0*dwdx49 + 1.0*dwdx50 - 1.0*dwdx51;
    JB[376] = -1.0*dwdx50;
    JB[379] = 1.0*dwdx51;
    JB[382] = 1.0*dwdx50;
    JB[384] = -1.0*dwdx52 + 1.0*dwdx54;
    JB[399] = 1.0*dwdx53;
    JB[400] = 1.0*dwdx52 - 1.0*dwdx53 - 1.0*dwdx54;
    JB[401] = 1.0*dwdx54;
    JB[406] = -1.0*dwdx52 + 1.0*dwdx53;
    JB[408] = 1.0*dwdx56;
    JB[424] = -1.0*dwdx56;
    JB[425] = 1.0*dwdx55 + 1.0*dwdx56 - 1.0*dwdx57;
    JB[427] = 1.0*dwdx57;
    JB[430] = -1.0*dwdx55 + 1.0*dwdx57;
    JB[449] = -1.0*dwdx60;
    JB[450] = 1.0*dwdx58 - 1.0*dwdx59 + 1.0*dwdx60;
    JB[453] = 1.0*dwdx59;
    JB[455] = -1.0*dwdx58 + 1.0*dwdx59;
    JB[456] = 1.0*dwdx62;
    JB[471] = -1.0*dwdx62;
    JB[473] = -1.0*dwdx63;
    JB[475] = 1.0*dwdx61 + 1.0*dwdx62 + 1.0*dwdx63 + 1.0*dwdx64;
    JB[477] = -1.0*dwdx64;
    JB[478] = 1.0*dwdx63;
    JB[499] = -1.0*dwdx66;
    JB[500] = 1.0*dwdx65;
    JB[522] = -1.0*dwdx68;
    JB[523] = 1.0*dwdx69;
    JB[525] = 1.0*dwdx67 + 1.0*dwdx68 - 1.0*dwdx69;
    JB[527] = 1.0*dwdx68;
    JB[529] = 1.0*dwdx70;
    JB[530] = -1.0*dwdx70;
    JB[531] = -1.0*dwdx73;
    JB[533] = 1.0*dwdx73;
    JB[536] = 1.0*dwdx71;
    JB[537] = -1.0*dwdx71;
    JB[538] = -1.0*dwdx74;
    JB[540] = 1.0*dwdx74;
    JB[543] = 1.0*dwdx72;
    JB[544] = -1.0*dwdx72;
    JB[545] = -1.0*dwdx75;
    JB[547] = 1.0*dwdx75;
    JB[550] = 1.0*dwdx70 + 1.0*dwdx71 + 1.0*dwdx72 + 1.0*dwdx73 + 1.0*dwdx74 + 1.0*dwdx75 + 1.0*dwdx76;
    JB[551] = -1.0*dwdx76;
    JB[556] = -1.0*dwdx77;
    JB[558] = -1.0*dwdx80;
    JB[559] = 1.0*dwdx77;
    JB[563] = -1.0*dwdx78;
    JB[566] = 1.0*dwdx78;
    JB[570] = -1.0*dwdx79;
    JB[573] = 1.0*dwdx79;
    JB[574] = 1.0*dwdx81;
    JB[575] = 1.0*dwdx77 + 1.0*dwdx78 + 1.0*dwdx79 - 1.0*dwdx81;
}