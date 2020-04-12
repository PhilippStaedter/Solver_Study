#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_ihekwaba1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3 - 1.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6;
    J[1] = 1.0*dwdx8 + 1.0*dwdx9;
    J[2] = 1.0*dwdx10 + 1.0*dwdx11;
    J[3] = 1.0*dwdx13 + 1.0*dwdx14;
    J[4] = 1.0*dwdx16 + 1.0*dwdx17;
    J[5] = 1.0*dwdx20 + 1.0*dwdx21;
    J[6] = 1.0*dwdx22 + 1.0*dwdx24;
    J[7] = -1.0*dwdx27;
    J[8] = -1.0*dwdx29;
    J[12] = -1.0*dwdx38;
    J[13] = -1.0*dwdx42;
    J[17] = -1.0*dwdx51;
    J[18] = -1.0*dwdx55;
    J[26] = 1.0*dwdx5;
    J[27] = -1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    J[28] = 1.0*dwdx12;
    J[33] = 1.0*dwdx27;
    J[48] = -1.0*dwdx64;
    J[52] = 1.0*dwdx2;
    J[53] = 1.0*dwdx7;
    J[54] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12;
    J[60] = 1.0*dwdx29;
    J[74] = 1.0*dwdx64;
    J[78] = 1.0*dwdx0;
    J[81] = -1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx15;
    J[82] = 1.0*dwdx18;
    J[90] = 1.0*dwdx38;
    J[100] = -1.0*dwdx70;
    J[104] = 1.0*dwdx3;
    J[107] = 1.0*dwdx15;
    J[108] = -1.0*dwdx16 - 1.0*dwdx17 - 1.0*dwdx18;
    J[117] = 1.0*dwdx42;
    J[126] = 1.0*dwdx70;
    J[130] = 1.0*dwdx1;
    J[135] = -1.0*dwdx19 - 1.0*dwdx20 - 1.0*dwdx21;
    J[136] = 1.0*dwdx23;
    J[147] = 1.0*dwdx51;
    J[152] = -1.0*dwdx65;
    J[156] = 1.0*dwdx4;
    J[161] = 1.0*dwdx19;
    J[162] = -1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24;
    J[174] = 1.0*dwdx55;
    J[178] = 1.0*dwdx65;
    J[182] = -1.0*dwdx5;
    J[183] = 1.0*dwdx9;
    J[189] = -1.0*dwdx25 - 1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28;
    J[190] = 1.0*dwdx30;
    J[191] = 1.0*dwdx33;
    J[193] = 1.0*dwdx37;
    J[204] = -1.0*dwdx66;
    J[208] = -1.0*dwdx2;
    J[210] = 1.0*dwdx11;
    J[215] = 1.0*dwdx25;
    J[216] = -1.0*dwdx29 - 1.0*dwdx30 - 1.0*dwdx31;
    J[218] = 1.0*dwdx35;
    J[230] = 1.0*dwdx66;
    J[241] = 1.0*dwdx26;
    J[243] = -1.0*dwdx32 - 1.0*dwdx33;
    J[244] = 1.0*dwdx34;
    J[257] = -1.0*dwdx72;
    J[269] = 1.0*dwdx32;
    J[270] = -1.0*dwdx34 - 1.0*dwdx35;
    J[283] = 1.0*dwdx72;
    J[297] = -1.0*dwdx36;
    J[309] = 1.0*dwdx71;
    J[311] = 1.0*dwdx76;
    J[312] = -1.0*dwdx0;
    J[315] = 1.0*dwdx14;
    J[324] = -1.0*dwdx38 - 1.0*dwdx39 - 1.0*dwdx40 - 1.0*dwdx41;
    J[325] = 1.0*dwdx43;
    J[326] = 1.0*dwdx46;
    J[328] = 1.0*dwdx50;
    J[334] = -1.0*dwdx67;
    J[338] = -1.0*dwdx3;
    J[342] = 1.0*dwdx16;
    J[350] = 1.0*dwdx39;
    J[351] = -1.0*dwdx42 - 1.0*dwdx43 - 1.0*dwdx44;
    J[353] = 1.0*dwdx48;
    J[360] = 1.0*dwdx67;
    J[376] = 1.0*dwdx40;
    J[378] = -1.0*dwdx45 - 1.0*dwdx46;
    J[379] = 1.0*dwdx47;
    J[387] = -1.0*dwdx73;
    J[404] = 1.0*dwdx45;
    J[405] = -1.0*dwdx47 - 1.0*dwdx48;
    J[413] = 1.0*dwdx73;
    J[432] = -1.0*dwdx49;
    J[441] = 1.0*dwdx77;
    J[442] = -1.0*dwdx1;
    J[447] = 1.0*dwdx21;
    J[459] = -1.0*dwdx51 - 1.0*dwdx52 - 1.0*dwdx53 - 1.0*dwdx54;
    J[460] = 1.0*dwdx56;
    J[461] = 1.0*dwdx59;
    J[463] = 1.0*dwdx63;
    J[464] = -1.0*dwdx68;
    J[468] = -1.0*dwdx4;
    J[474] = 1.0*dwdx24;
    J[485] = 1.0*dwdx52;
    J[486] = -1.0*dwdx55 - 1.0*dwdx56 - 1.0*dwdx57;
    J[488] = 1.0*dwdx61;
    J[490] = 1.0*dwdx68;
    J[511] = 1.0*dwdx53;
    J[513] = -1.0*dwdx58 - 1.0*dwdx59;
    J[514] = 1.0*dwdx60;
    J[517] = -1.0*dwdx74;
    J[539] = 1.0*dwdx58;
    J[540] = -1.0*dwdx60 - 1.0*dwdx61;
    J[543] = 1.0*dwdx74;
    J[567] = -1.0*dwdx62;
    J[571] = 1.0*dwdx78;
    J[573] = -1.0*dwdx7;
    J[574] = 1.0*dwdx10 + 1.0*dwdx12;
    J[575] = -1.0*dwdx15;
    J[576] = 1.0*dwdx17 + 1.0*dwdx18;
    J[577] = -1.0*dwdx19;
    J[578] = 1.0*dwdx22 + 1.0*dwdx23;
    J[579] = -1.0*dwdx25;
    J[580] = 1.0*dwdx30 + 1.0*dwdx31;
    J[584] = -1.0*dwdx39;
    J[585] = 1.0*dwdx43 + 1.0*dwdx44;
    J[589] = -1.0*dwdx52;
    J[590] = 1.0*dwdx56 + 1.0*dwdx57;
    J[594] = -1.0*dwdx64 - 1.0*dwdx65 - 1.0*dwdx66 - 1.0*dwdx67 - 1.0*dwdx68 - 1.0*dwdx69 - 1.0*dwdx70;
    J[595] = 1.0*dwdx75;
    J[607] = -1.0*dwdx32;
    J[608] = 1.0*dwdx34;
    J[612] = -1.0*dwdx45;
    J[613] = 1.0*dwdx47;
    J[617] = -1.0*dwdx58;
    J[618] = 1.0*dwdx60;
    J[620] = 1.0*dwdx69;
    J[621] = -1.0*dwdx72 - 1.0*dwdx73 - 1.0*dwdx74 - 1.0*dwdx75;
}