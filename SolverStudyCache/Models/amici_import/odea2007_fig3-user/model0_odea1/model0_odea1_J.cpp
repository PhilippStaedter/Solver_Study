#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_odea1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3 - 1.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6;
    J[1] = 1.0*dwdx7 - 1.0*dwdx9;
    J[2] = 1.0*dwdx10 - 1.0*dwdx12;
    J[3] = -1.0*dwdx14;
    J[5] = -1.0*dwdx20;
    J[8] = 1.0*dwdx28 - 1.0*dwdx30;
    J[9] = 1.0*dwdx31 - 1.0*dwdx33;
    J[10] = -1.0*dwdx35;
    J[12] = -1.0*dwdx41;
    J[15] = 1.0*dwdx49 - 1.0*dwdx51;
    J[16] = 1.0*dwdx52 - 1.0*dwdx54;
    J[17] = -1.0*dwdx56;
    J[19] = -1.0*dwdx62;
    J[24] = 1.0*dwdx4;
    J[25] = -1.0*dwdx7 - 1.0*dwdx8 + 1.0*dwdx9;
    J[26] = -1.0*dwdx11;
    J[29] = 1.0*dwdx20;
    J[46] = -1.0*dwdx70;
    J[48] = 1.0*dwdx1;
    J[49] = 1.0*dwdx8;
    J[50] = -1.0*dwdx10 + 1.0*dwdx11 + 1.0*dwdx12;
    J[51] = 1.0*dwdx14;
    J[70] = 1.0*dwdx70;
    J[72] = -1.0*dwdx1;
    J[74] = -1.0*dwdx12;
    J[75] = -1.0*dwdx13 - 1.0*dwdx14 + 1.0*dwdx15;
    J[76] = 1.0*dwdx18;
    J[77] = 1.0*dwdx21;
    J[94] = 1.0*dwdx73;
    J[100] = -1.0*dwdx16 + 1.0*dwdx17 - 1.0*dwdx18;
    J[103] = 1.0*dwdx26;
    J[119] = 1.0*dwdx77;
    J[120] = -1.0*dwdx4;
    J[121] = -1.0*dwdx9;
    J[123] = -1.0*dwdx15;
    J[125] = -1.0*dwdx19 - 1.0*dwdx20 - 1.0*dwdx21 - 1.0*dwdx22;
    J[126] = 1.0*dwdx24;
    J[127] = -1.0*dwdx27;
    J[142] = -1.0*dwdx73;
    J[150] = -1.0*dwdx23;
    J[167] = 1.0*dwdx80;
    J[172] = -1.0*dwdx17;
    J[173] = 1.0*dwdx22;
    J[175] = -1.0*dwdx25 - 1.0*dwdx26 + 1.0*dwdx27;
    J[191] = -1.0*dwdx77;
    J[192] = 1.0*dwdx5;
    J[200] = -1.0*dwdx28 - 1.0*dwdx29 + 1.0*dwdx30;
    J[201] = -1.0*dwdx32;
    J[204] = 1.0*dwdx41;
    J[214] = -1.0*dwdx71;
    J[216] = 1.0*dwdx2;
    J[224] = 1.0*dwdx29;
    J[225] = -1.0*dwdx31 + 1.0*dwdx32 + 1.0*dwdx33;
    J[226] = 1.0*dwdx35;
    J[238] = 1.0*dwdx71;
    J[240] = -1.0*dwdx2;
    J[249] = -1.0*dwdx33;
    J[250] = -1.0*dwdx34 - 1.0*dwdx35 + 1.0*dwdx36;
    J[251] = 1.0*dwdx39;
    J[252] = 1.0*dwdx42;
    J[262] = 1.0*dwdx74;
    J[275] = -1.0*dwdx37 + 1.0*dwdx38 - 1.0*dwdx39;
    J[278] = 1.0*dwdx47;
    J[287] = 1.0*dwdx78;
    J[288] = -1.0*dwdx5;
    J[296] = -1.0*dwdx30;
    J[298] = -1.0*dwdx36;
    J[300] = -1.0*dwdx40 - 1.0*dwdx41 - 1.0*dwdx42 - 1.0*dwdx43;
    J[301] = 1.0*dwdx45;
    J[302] = -1.0*dwdx48;
    J[310] = -1.0*dwdx74;
    J[325] = -1.0*dwdx44;
    J[347] = -1.0*dwdx38;
    J[348] = 1.0*dwdx43;
    J[350] = -1.0*dwdx46 - 1.0*dwdx47 + 1.0*dwdx48;
    J[359] = -1.0*dwdx78;
    J[360] = 1.0*dwdx6;
    J[375] = -1.0*dwdx49 - 1.0*dwdx50 + 1.0*dwdx51;
    J[376] = -1.0*dwdx53;
    J[379] = 1.0*dwdx62;
    J[382] = -1.0*dwdx72;
    J[384] = 1.0*dwdx3;
    J[399] = 1.0*dwdx50;
    J[400] = -1.0*dwdx52 + 1.0*dwdx53 + 1.0*dwdx54;
    J[401] = 1.0*dwdx56;
    J[406] = 1.0*dwdx72;
    J[408] = -1.0*dwdx3;
    J[424] = -1.0*dwdx54;
    J[425] = -1.0*dwdx55 - 1.0*dwdx56 + 1.0*dwdx57;
    J[426] = 1.0*dwdx60;
    J[427] = 1.0*dwdx63;
    J[430] = 1.0*dwdx75;
    J[450] = -1.0*dwdx58 + 1.0*dwdx59 - 1.0*dwdx60;
    J[453] = 1.0*dwdx68;
    J[455] = 1.0*dwdx79;
    J[456] = -1.0*dwdx6;
    J[471] = -1.0*dwdx51;
    J[473] = -1.0*dwdx57;
    J[475] = -1.0*dwdx61 - 1.0*dwdx62 - 1.0*dwdx63 - 1.0*dwdx64;
    J[476] = 1.0*dwdx66;
    J[477] = -1.0*dwdx69;
    J[478] = -1.0*dwdx75;
    J[500] = -1.0*dwdx65;
    J[522] = -1.0*dwdx59;
    J[523] = 1.0*dwdx64;
    J[525] = -1.0*dwdx67 - 1.0*dwdx68 + 1.0*dwdx69;
    J[527] = -1.0*dwdx79;
    J[529] = -1.0*dwdx8;
    J[530] = 1.0*dwdx10 - 1.0*dwdx11;
    J[531] = 1.0*dwdx13 - 1.0*dwdx15;
    J[533] = -1.0*dwdx21;
    J[536] = -1.0*dwdx29;
    J[537] = 1.0*dwdx31 - 1.0*dwdx32;
    J[538] = 1.0*dwdx34 - 1.0*dwdx36;
    J[540] = -1.0*dwdx42;
    J[543] = -1.0*dwdx50;
    J[544] = 1.0*dwdx52 - 1.0*dwdx53;
    J[545] = 1.0*dwdx55 - 1.0*dwdx57;
    J[547] = -1.0*dwdx63;
    J[550] = -1.0*dwdx70 - 1.0*dwdx71 - 1.0*dwdx72 - 1.0*dwdx73 - 1.0*dwdx74 - 1.0*dwdx75 - 1.0*dwdx76;
    J[551] = -1.0*dwdx81;
    J[556] = 1.0*dwdx16 - 1.0*dwdx17;
    J[559] = -1.0*dwdx26;
    J[563] = 1.0*dwdx37 - 1.0*dwdx38;
    J[566] = -1.0*dwdx47;
    J[570] = 1.0*dwdx58 - 1.0*dwdx59;
    J[573] = -1.0*dwdx68;
    J[574] = 1.0*dwdx76;
    J[575] = -1.0*dwdx77 - 1.0*dwdx78 - 1.0*dwdx79 + 1.0*dwdx81;
}