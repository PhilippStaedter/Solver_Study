#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_Sengupta2015(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0;
    JDiag[2] = -1.0*dwdx1 - 1.0*dwdx2;
    JDiag[3] = -1.0*dwdx3;
    JDiag[6] = -1.0*dwdx4 - 1.0*dwdx5;
    JDiag[12] = -1.0*dwdx6;
    JDiag[18] = -1.0*dwdx7;
    JDiag[21] = -1.0*dwdx8;
    JDiag[24] = -1.0*dwdx9;
    JDiag[25] = -1.0*dwdx10;
    JDiag[26] = -1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx13 - 1.0*dwdx14;
    JDiag[27] = -1.0*dwdx15;
    JDiag[28] = -1.0*dwdx16 - 1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19;
    JDiag[29] = -1.0*dwdx20;
    JDiag[43] = -1.0*dwdx21;
    JDiag[44] = -1.0*dwdx22;
    JDiag[45] = -1.0*dwdx23;
    JDiag[46] = -1.0*dwdx24;
    JDiag[50] = -1.0*dwdx25;
    JDiag[51] = -1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28 - 1.0*dwdx29;
    JDiag[52] = -1.0*dwdx30;
    JDiag[53] = -1.0*dwdx31;
    JDiag[58] = -1.0*dwdx36;
    JDiag[59] = -1.0*dwdx37;
    JDiag[60] = -1.0*dwdx38;
    JDiag[62] = -1.0*dwdx39;
    JDiag[65] = -1.0*dwdx40;
    JDiag[67] = -1.0*dwdx41;
    JDiag[69] = -1.0*dwdx42;
    JDiag[71] = -1.0*dwdx43 - 1.0*dwdx44;
    JDiag[80] = -1.0*dwdx45;
    JDiag[81] = -1.0*dwdx46;
    JDiag[82] = -1.0*dwdx47;
    JDiag[83] = -1.0*dwdx48 - 1.0*dwdx49;
    JDiag[84] = -1.0*dwdx50 - 1.0*dwdx51;
    JDiag[85] = -1.0*dwdx52;
    JDiag[87] = -1.0*dwdx53;
    JDiag[88] = -1.0*dwdx54;
    JDiag[89] = -1.0*dwdx55;
    JDiag[90] = -1.0*dwdx56;
    JDiag[91] = -1.0*dwdx57;
    JDiag[96] = -1.0*dwdx58 - 1.0*dwdx59 - 1.0*dwdx60;
    JDiag[98] = -1.0*dwdx61 - 1.0*dwdx62 - 1.0*dwdx63 - 1.0*dwdx64;
    JDiag[100] = -1.0*dwdx65;
    JDiag[102] = -1.0*dwdx66;
    JDiag[104] = -1.0*dwdx67;
    JDiag[106] = -1.0*dwdx68;
    JDiag[107] = -1.0*dwdx69;
    JDiag[108] = -1.0*dwdx71 - 1.0*dwdx72 - 1.0*dwdx73;
    JDiag[109] = -1.0*dwdx74;
    JDiag[111] = -1.0*dwdx77;
    JDiag[112] = -1.0*dwdx78;
    JDiag[113] = -1.0*dwdx79;
    JDiag[114] = -1.0*dwdx80 - 1.0*dwdx81 - 1.0*dwdx82;
    JDiag[121] = -1.0*dwdx83;
    JDiag[124] = -1.0*dwdx84;
    JDiag[126] = -1.0*dwdx85;
    JDiag[127] = -1.0*dwdx86;
    JDiag[131] = -1.0*dwdx87;
    JDiag[132] = -1.0*dwdx88;
    JDiag[134] = -1.0*dwdx90;
    JDiag[135] = -1.0*dwdx91 - 1.0*dwdx92 - 1.0*dwdx93 - 1.0*dwdx94 - 1.0*dwdx95;
    JDiag[137] = -1.0*dwdx96;
    JDiag[141] = -1.0*dwdx97;
    JDiag[142] = -1.0*dwdx98;
    JDiag[143] = -1.0*dwdx100 - 1.0*dwdx99;
    JDiag[144] = -1.0*dwdx101;
    JDiag[145] = -1.0*dwdx102;
    JDiag[146] = -1.0*dwdx103;
    JDiag[147] = -1.0*dwdx104 - 1.0*dwdx105;
    JDiag[148] = -1.0*dwdx106;
    JDiag[149] = -1.0*dwdx107;
    JDiag[150] = -1.0*dwdx108;
    JDiag[151] = -1.0*dwdx109;
    JDiag[152] = -1.0*dwdx110;
    JDiag[153] = -1.0*dwdx111;
    JDiag[154] = -1.0*dwdx112;
    JDiag[155] = -1.0*dwdx113;
    JDiag[156] = -1.0*dwdx114;
    JDiag[157] = -1.0*dwdx115 - 1.0*dwdx116;
    JDiag[158] = -1.0*dwdx117;
    JDiag[159] = -1.0*dwdx118;
    JDiag[160] = -1.0*dwdx119;
    JDiag[161] = -1.0*dwdx120;
    JDiag[162] = -1.0*dwdx121;
    JDiag[163] = -1.0*dwdx122;
    JDiag[164] = -1.0*dwdx123 - 1.0*dwdx124;
    JDiag[165] = -1.0*dwdx125;
    JDiag[166] = -1.0*dwdx126;
    JDiag[167] = -1.0*dwdx127;
    JDiag[168] = -1.0*dwdx128;
    JDiag[169] = -1.0*dwdx129;
    JDiag[170] = -1.0*dwdx130;
    JDiag[171] = -1.0*dwdx131 - 1.0*dwdx132;
    JDiag[172] = -1.0*dwdx133;
    JDiag[173] = -1.0*dwdx134;
    JDiag[174] = -1.0*dwdx135;
    JDiag[175] = -1.0*dwdx136;
    JDiag[176] = -1.0*dwdx137;
    JDiag[177] = -1.0*dwdx138;
    JDiag[178] = -1.0*dwdx139 - 1.0*dwdx140;
    JDiag[179] = -1.0*dwdx141;
    JDiag[180] = -1.0*dwdx142;
    JDiag[181] = -1.0*dwdx143;
    JDiag[182] = -1.0*dwdx144;
    JDiag[183] = -1.0*dwdx145;
    JDiag[184] = -1.0*dwdx146;
    JDiag[185] = -1.0*dwdx147 - 1.0*dwdx148;
    JDiag[186] = -1.0*dwdx149;
    JDiag[187] = -1.0*dwdx150;
    JDiag[188] = -1.0*dwdx151;
    JDiag[189] = -1.0*dwdx152;
    JDiag[190] = -1.0*dwdx153 - 1.0*dwdx154;
    JDiag[191] = -1.0*dwdx155;
    JDiag[192] = -1.0*dwdx156;
    JDiag[193] = -1.0*dwdx157;
    JDiag[194] = -1.0*dwdx158 - 1.0*dwdx159;
    JDiag[195] = -1.0*dwdx160;
    JDiag[199] = -1.0*dwdx161;
    JDiag[200] = -1.0*dwdx162;
    JDiag[201] = -1.0*dwdx163;
    JDiag[202] = -1.0*dwdx165;
    JDiag[205] = -1.0*dwdx166;
    JDiag[206] = -1.0*dwdx167;
    JDiag[207] = -1.0*dwdx168;
    JDiag[209] = -1.0*dwdx169;
    JDiag[214] = -1.0*dwdx170 - 1.0*dwdx171;
    JDiag[216] = -1.0*dwdx172 - 1.0*dwdx173 - 1.0*dwdx174;
    JDiag[217] = -1.0*dwdx175 - 1.0*dwdx176;
    JDiag[218] = -1.0*dwdx177 - 1.0*dwdx178;
    JDiag[222] = -1.0*dwdx179 - 1.0*dwdx180;
    JDiag[224] = -1.0*dwdx181;
    JDiag[225] = -1.0*dwdx182;
    JDiag[226] = -1.0*dwdx183 - 1.0*dwdx184;
    JDiag[227] = -1.0*dwdx185;
    JDiag[230] = -1.0*dwdx186 - 1.0*dwdx187 - 1.0*dwdx188;
    JDiag[234] = -1.0*dwdx189;
    JDiag[235] = -1.0*dwdx190;
    JDiag[236] = -1.0*dwdx191 - 1.0*dwdx192;
    JDiag[237] = -1.0*dwdx193;
    JDiag[238] = -1.0*dwdx194;
}