#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_ihekwaba1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0 + 1.0*dwdx1 + 1.0*dwdx2 + 1.0*dwdx3 + 1.0*dwdx4 + 1.0*dwdx5 + 1.0*dwdx6;
    JB[1] = -1.0*dwdx5;
    JB[2] = -1.0*dwdx2;
    JB[3] = -1.0*dwdx0;
    JB[4] = -1.0*dwdx3;
    JB[5] = -1.0*dwdx1;
    JB[6] = -1.0*dwdx4;
    JB[7] = 1.0*dwdx5;
    JB[8] = 1.0*dwdx2;
    JB[12] = 1.0*dwdx0;
    JB[13] = 1.0*dwdx3;
    JB[17] = 1.0*dwdx1;
    JB[18] = 1.0*dwdx4;
    JB[26] = -1.0*dwdx8 - 1.0*dwdx9;
    JB[27] = 1.0*dwdx7 + 1.0*dwdx8 + 1.0*dwdx9;
    JB[28] = -1.0*dwdx7;
    JB[33] = -1.0*dwdx9;
    JB[48] = 1.0*dwdx7;
    JB[52] = -1.0*dwdx10 - 1.0*dwdx11;
    JB[53] = -1.0*dwdx12;
    JB[54] = 1.0*dwdx10 + 1.0*dwdx11 + 1.0*dwdx12;
    JB[60] = -1.0*dwdx11;
    JB[74] = -1.0*dwdx10 - 1.0*dwdx12;
    JB[78] = -1.0*dwdx13 - 1.0*dwdx14;
    JB[81] = 1.0*dwdx13 + 1.0*dwdx14 + 1.0*dwdx15;
    JB[82] = -1.0*dwdx15;
    JB[90] = -1.0*dwdx14;
    JB[100] = 1.0*dwdx15;
    JB[104] = -1.0*dwdx16 - 1.0*dwdx17;
    JB[107] = -1.0*dwdx18;
    JB[108] = 1.0*dwdx16 + 1.0*dwdx17 + 1.0*dwdx18;
    JB[117] = -1.0*dwdx16;
    JB[126] = -1.0*dwdx17 - 1.0*dwdx18;
    JB[130] = -1.0*dwdx20 - 1.0*dwdx21;
    JB[135] = 1.0*dwdx19 + 1.0*dwdx20 + 1.0*dwdx21;
    JB[136] = -1.0*dwdx19;
    JB[147] = -1.0*dwdx21;
    JB[152] = 1.0*dwdx19;
    JB[156] = -1.0*dwdx22 - 1.0*dwdx24;
    JB[161] = -1.0*dwdx23;
    JB[162] = 1.0*dwdx22 + 1.0*dwdx23 + 1.0*dwdx24;
    JB[174] = -1.0*dwdx24;
    JB[178] = -1.0*dwdx22 - 1.0*dwdx23;
    JB[182] = 1.0*dwdx27;
    JB[183] = -1.0*dwdx27;
    JB[189] = 1.0*dwdx25 + 1.0*dwdx26 + 1.0*dwdx27 + 1.0*dwdx28;
    JB[190] = -1.0*dwdx25;
    JB[191] = -1.0*dwdx26;
    JB[204] = 1.0*dwdx25;
    JB[208] = 1.0*dwdx29;
    JB[210] = -1.0*dwdx29;
    JB[215] = -1.0*dwdx30;
    JB[216] = 1.0*dwdx29 + 1.0*dwdx30 + 1.0*dwdx31;
    JB[230] = -1.0*dwdx30 - 1.0*dwdx31;
    JB[241] = -1.0*dwdx33;
    JB[243] = 1.0*dwdx32 + 1.0*dwdx33;
    JB[244] = -1.0*dwdx32;
    JB[257] = 1.0*dwdx32;
    JB[268] = -1.0*dwdx35;
    JB[269] = -1.0*dwdx34;
    JB[270] = 1.0*dwdx34 + 1.0*dwdx35;
    JB[283] = -1.0*dwdx34;
    JB[293] = -1.0*dwdx37;
    JB[297] = 1.0*dwdx36;
    JB[312] = 1.0*dwdx38;
    JB[315] = -1.0*dwdx38;
    JB[324] = 1.0*dwdx38 + 1.0*dwdx39 + 1.0*dwdx40 + 1.0*dwdx41;
    JB[325] = -1.0*dwdx39;
    JB[326] = -1.0*dwdx40;
    JB[334] = 1.0*dwdx39;
    JB[338] = 1.0*dwdx42;
    JB[342] = -1.0*dwdx42;
    JB[350] = -1.0*dwdx43;
    JB[351] = 1.0*dwdx42 + 1.0*dwdx43 + 1.0*dwdx44;
    JB[360] = -1.0*dwdx43 - 1.0*dwdx44;
    JB[376] = -1.0*dwdx46;
    JB[378] = 1.0*dwdx45 + 1.0*dwdx46;
    JB[379] = -1.0*dwdx45;
    JB[387] = 1.0*dwdx45;
    JB[403] = -1.0*dwdx48;
    JB[404] = -1.0*dwdx47;
    JB[405] = 1.0*dwdx47 + 1.0*dwdx48;
    JB[413] = -1.0*dwdx47;
    JB[428] = -1.0*dwdx50;
    JB[432] = 1.0*dwdx49;
    JB[442] = 1.0*dwdx51;
    JB[447] = -1.0*dwdx51;
    JB[459] = 1.0*dwdx51 + 1.0*dwdx52 + 1.0*dwdx53 + 1.0*dwdx54;
    JB[460] = -1.0*dwdx52;
    JB[461] = -1.0*dwdx53;
    JB[464] = 1.0*dwdx52;
    JB[468] = 1.0*dwdx55;
    JB[474] = -1.0*dwdx55;
    JB[485] = -1.0*dwdx56;
    JB[486] = 1.0*dwdx55 + 1.0*dwdx56 + 1.0*dwdx57;
    JB[490] = -1.0*dwdx56 - 1.0*dwdx57;
    JB[511] = -1.0*dwdx59;
    JB[513] = 1.0*dwdx58 + 1.0*dwdx59;
    JB[514] = -1.0*dwdx58;
    JB[517] = 1.0*dwdx58;
    JB[538] = -1.0*dwdx61;
    JB[539] = -1.0*dwdx60;
    JB[540] = 1.0*dwdx60 + 1.0*dwdx61;
    JB[543] = -1.0*dwdx60;
    JB[563] = -1.0*dwdx63;
    JB[567] = 1.0*dwdx62;
    JB[573] = 1.0*dwdx64;
    JB[574] = -1.0*dwdx64;
    JB[575] = 1.0*dwdx70;
    JB[576] = -1.0*dwdx70;
    JB[577] = 1.0*dwdx65;
    JB[578] = -1.0*dwdx65;
    JB[579] = 1.0*dwdx66;
    JB[580] = -1.0*dwdx66;
    JB[584] = 1.0*dwdx67;
    JB[585] = -1.0*dwdx67;
    JB[589] = 1.0*dwdx68;
    JB[590] = -1.0*dwdx68;
    JB[594] = 1.0*dwdx64 + 1.0*dwdx65 + 1.0*dwdx66 + 1.0*dwdx67 + 1.0*dwdx68 + 1.0*dwdx69 + 1.0*dwdx70;
    JB[595] = -1.0*dwdx69;
    JB[607] = 1.0*dwdx72;
    JB[608] = -1.0*dwdx72;
    JB[609] = -1.0*dwdx71;
    JB[612] = 1.0*dwdx73;
    JB[613] = -1.0*dwdx73;
    JB[617] = 1.0*dwdx74;
    JB[618] = -1.0*dwdx74;
    JB[620] = -1.0*dwdx75;
    JB[621] = 1.0*dwdx72 + 1.0*dwdx73 + 1.0*dwdx74 + 1.0*dwdx75;
    JB[661] = -1.0*dwdx76;
    JB[666] = -1.0*dwdx77;
    JB[671] = -1.0*dwdx78;
}