#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Qi2013a(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[3] = 1.0*Source;
            break;
        case 1:
            dwdp[4] = 1.0*NatP*ROS;
            break;
        case 2:
            dwdp[5] = 1.0*MisP;
            break;
        case 3:
            dwdp[6] = 1.0*E3*MisP;
            break;
        case 4:
            dwdp[7] = 1.0*E3_MisP;
            break;
        case 5:
            dwdp[8] = 1.0*ATP*E1*Ub/(1.0*ATP + 5000);
            break;
        case 6:
            dwdp[9] = 1.0*E1_Ub*E2;
            break;
        case 7:
            dwdp[10] = 1.0*E2_Ub*E3_MisP;
            dwdp[105] = 1.0*E2_Ub*E3SUB_SUB_misfolded;
            dwdp[152] = 1.0*E2_Ub*Parkin_asyn_dam;
            break;
        case 8:
            dwdp[11] = 1.0*E2_Ub*E3_MisP_Ub;
            dwdp[12] = 1.0*E2_Ub*E3_MisP_Ub2;
            dwdp[13] = 1.0*E2_Ub*E3_MisP_Ub3;
            dwdp[14] = 1.0*E2_Ub*E3_MisP_Ub4;
            dwdp[15] = 1.0*E2_Ub*E3_MisP_Ub5;
            dwdp[16] = 1.0*E2_Ub*E3_MisP_Ub6;
            dwdp[17] = 1.0*E2_Ub*E3_MisP_Ub7;
            dwdp[106] = 1.0*E2_Ub*E3SUB_SUB_misfolded_Ub;
            dwdp[107] = 1.0*E2_Ub*E3SUB_SUB_misfolded_Ub2;
            dwdp[108] = 1.0*E2_Ub*E3SUB_SUB_misfolded_Ub3;
            dwdp[109] = 1.0*E2_Ub*E3SUB_SUB_misfolded_Ub4;
            dwdp[110] = 1.0*E2_Ub*E3SUB_SUB_misfolded_Ub5;
            dwdp[111] = 1.0*E2_Ub*E3SUB_SUB_misfolded_Ub6;
            dwdp[112] = 1.0*E2_Ub*E3SUB_SUB_misfolded_Ub7;
            dwdp[153] = 1.0*E2_Ub*Parkin_asyn_dam_Ub;
            dwdp[154] = 1.0*E2_Ub*Parkin_asyn_dam_Ub2;
            dwdp[155] = 1.0*E2_Ub*Parkin_asyn_dam_Ub3;
            dwdp[156] = 1.0*E2_Ub*Parkin_asyn_dam_Ub4;
            dwdp[157] = 1.0*E2_Ub*Parkin_asyn_dam_Ub5;
            dwdp[158] = 1.0*E2_Ub*Parkin_asyn_dam_Ub6;
            dwdp[159] = 1.0*E2_Ub*Parkin_asyn_dam_Ub7;
            break;
        case 9:
            dwdp[26] = 1.0*E3_MisP_Ub8_DUB;
            dwdp[27] = 1.0*E3_MisP_Ub7_DUB;
            dwdp[28] = 1.0*E3_MisP_Ub6_DUB;
            dwdp[29] = 1.0*E3_MisP_Ub5_DUB;
            dwdp[30] = 1.0*E3_MisP_Ub4_DUB;
            dwdp[31] = 1.0*E3_MisP_Ub3_DUB;
            dwdp[32] = 1.0*E3_MisP_Ub2_DUB;
            dwdp[33] = 1.0*E3_MisP_Ub_DUB;
            dwdp[168] = 1.0*Parkin_asyn_dam_Ub8_DUB;
            dwdp[169] = 1.0*Parkin_asyn_dam_Ub7_DUB;
            dwdp[170] = 1.0*Parkin_asyn_dam_Ub6_DUB;
            dwdp[171] = 1.0*Parkin_asyn_dam_Ub5_DUB;
            dwdp[172] = 1.0*Parkin_asyn_dam_Ub4_DUB;
            dwdp[173] = 1.0*Parkin_asyn_dam_Ub3_DUB;
            dwdp[174] = 1.0*Parkin_asyn_dam_Ub2_DUB;
            dwdp[175] = 1.0*Parkin_asyn_dam_Ub_DUB;
            break;
        case 10:
            dwdp[34] = 1.0*E3_MisP_Ub4*Proteasome;
            dwdp[35] = 1.0*E3_MisP_Ub5*Proteasome;
            dwdp[36] = 1.0*E3_MisP_Ub6*Proteasome;
            dwdp[37] = 1.0*E3_MisP_Ub7*Proteasome;
            dwdp[38] = 1.0*E3_MisP_Ub8*Proteasome;
            dwdp[129] = 1.0*E3SUB_SUB_misfolded_Ub4*Proteasome;
            dwdp[130] = 1.0*E3SUB_SUB_misfolded_Ub5*Proteasome;
            dwdp[131] = 1.0*E3SUB_SUB_misfolded_Ub6*Proteasome;
            dwdp[132] = 1.0*E3SUB_SUB_misfolded_Ub7*Proteasome;
            dwdp[133] = 1.0*E3SUB_SUB_misfolded_Ub8*Proteasome;
            dwdp[176] = 1.0*Parkin_asyn_dam_Ub4*Proteasome;
            dwdp[177] = 1.0*Parkin_asyn_dam_Ub5*Proteasome;
            dwdp[178] = 1.0*Parkin_asyn_dam_Ub6*Proteasome;
            dwdp[179] = 1.0*Parkin_asyn_dam_Ub7*Proteasome;
            dwdp[180] = 1.0*Parkin_asyn_dam_Ub8*Proteasome;
            break;
        case 11:
            dwdp[39] = 1.0*DUB*MisP_Ub8_Proteasome;
            dwdp[40] = 1.0*DUB*MisP_Ub7_Proteasome;
            dwdp[41] = 1.0*DUB*MisP_Ub6_Proteasome;
            dwdp[42] = 1.0*DUB*MisP_Ub5_Proteasome;
            dwdp[43] = 1.0*DUB*MisP_Ub4_Proteasome;
            dwdp[134] = 1.0*DUB*SUB_misfolded_Ub8_Proteasome;
            dwdp[135] = 1.0*DUB*SUB_misfolded_Ub7_Proteasome;
            dwdp[136] = 1.0*DUB*SUB_misfolded_Ub6_Proteasome;
            dwdp[137] = 1.0*DUB*SUB_misfolded_Ub5_Proteasome;
            dwdp[138] = 1.0*DUB*SUB_misfolded_Ub4_Proteasome;
            dwdp[181] = 1.0*DUB*asyn_dam_Ub8_Proteasome;
            dwdp[182] = 1.0*DUB*asyn_dam_Ub7_Proteasome;
            dwdp[183] = 1.0*DUB*asyn_dam_Ub6_Proteasome;
            dwdp[184] = 1.0*DUB*asyn_dam_Ub5_Proteasome;
            dwdp[185] = 1.0*DUB*asyn_dam_Ub4_Proteasome;
            break;
        case 12:
            dwdp[44] = 1.0*ATP*MisP_Ub4_Proteasome*kproteff/(1.0*ATP + 5000);
            dwdp[45] = 1.0*ATP*MisP_Ub5_Proteasome*kproteff/(1.0*ATP + 5000);
            dwdp[46] = 1.0*ATP*MisP_Ub6_Proteasome*kproteff/(1.0*ATP + 5000);
            dwdp[47] = 1.0*ATP*MisP_Ub7_Proteasome*kproteff/(1.0*ATP + 5000);
            dwdp[48] = 1.0*ATP*MisP_Ub8_Proteasome*kproteff/(1.0*ATP + 5000);
            dwdp[139] = 1.0*ATP*SUB_misfolded_Ub4_Proteasome*kproteff/(1.0*ATP + 5000);
            dwdp[140] = 1.0*ATP*SUB_misfolded_Ub5_Proteasome*kproteff/(1.0*ATP + 5000);
            dwdp[141] = 1.0*ATP*SUB_misfolded_Ub6_Proteasome*kproteff/(1.0*ATP + 5000);
            dwdp[142] = 1.0*ATP*SUB_misfolded_Ub7_Proteasome*kproteff/(1.0*ATP + 5000);
            dwdp[143] = 1.0*ATP*SUB_misfolded_Ub8_Proteasome*kproteff/(1.0*ATP + 5000);
            dwdp[186] = 1.0*ATP*asyn_dam_Ub4_Proteasome*kproteff/(1.0*ATP + 5000);
            dwdp[187] = 1.0*ATP*asyn_dam_Ub5_Proteasome*kproteff/(1.0*ATP + 5000);
            dwdp[188] = 1.0*ATP*asyn_dam_Ub6_Proteasome*kproteff/(1.0*ATP + 5000);
            dwdp[189] = 1.0*ATP*asyn_dam_Ub7_Proteasome*kproteff/(1.0*ATP + 5000);
            dwdp[190] = 1.0*ATP*asyn_dam_Ub8_Proteasome*kproteff/(1.0*ATP + 5000);
            break;
        case 13:
            dwdp[49] = 0.5*MisP*(1.0*MisP - 1);
            break;
        case 14:
            dwdp[50] = 1.0*AggP1*MisP;
            dwdp[51] = 1.0*AggP2*MisP;
            dwdp[52] = 1.0*AggP3*MisP;
            dwdp[53] = 1.0*AggP4*MisP;
            dwdp[59] = 1.0*AggP5*MisP;
            break;
        case 15:
            dwdp[58] = 1.0*AggP1;
            break;
        case 16:
            dwdp[57] = 1.0*AggP2;
            break;
        case 17:
            dwdp[56] = 1.0*AggP3;
            break;
        case 18:
            dwdp[55] = 1.0*AggP4;
            break;
        case 19:
            dwdp[54] = 1.0*AggP5;
            break;
        case 20:
            dwdp[60] = 1.0*MisP*SeqAggP;
            dwdp[207] = 1.0*SeqAggP*asyn;
            dwdp[229] = 1.0*SeqAggP*asyn_dam;
            dwdp[268] = 1.0*SeqAggP*UCHL1_damaged;
            dwdp[269] = 1.0*Lamp2a_UCHL1_damaged*SeqAggP;
            dwdp[291] = 1.0*SUB_misfolded*SeqAggP;
            break;
        case 21:
            dwdp[61] = 1.0*E3_MisP*SeqAggP;
            dwdp[62] = 1.0*E3_MisP_Ub*SeqAggP;
            dwdp[63] = 1.0*E3_MisP_Ub2*SeqAggP;
            dwdp[64] = 1.0*E3_MisP_Ub3*SeqAggP;
            dwdp[65] = 1.0*E3_MisP_Ub4*SeqAggP;
            dwdp[66] = 1.0*E3_MisP_Ub5*SeqAggP;
            dwdp[67] = 1.0*E3_MisP_Ub6*SeqAggP;
            dwdp[68] = 1.0*E3_MisP_Ub7*SeqAggP;
            dwdp[69] = 1.0*E3_MisP_Ub8*SeqAggP;
            dwdp[70] = 1.0*E3_MisP_Ub_DUB*SeqAggP;
            dwdp[71] = 1.0*E3_MisP_Ub2_DUB*SeqAggP;
            dwdp[72] = 1.0*E3_MisP_Ub3_DUB*SeqAggP;
            dwdp[73] = 1.0*E3_MisP_Ub4_DUB*SeqAggP;
            dwdp[74] = 1.0*E3_MisP_Ub5_DUB*SeqAggP;
            dwdp[75] = 1.0*E3_MisP_Ub6_DUB*SeqAggP;
            dwdp[76] = 1.0*E3_MisP_Ub7_DUB*SeqAggP;
            dwdp[77] = 1.0*E3_MisP_Ub8_DUB*SeqAggP;
            dwdp[230] = 1.0*Parkin_asyn_dam*SeqAggP;
            dwdp[231] = 1.0*Parkin_asyn_dam_Ub*SeqAggP;
            dwdp[232] = 1.0*Parkin_asyn_dam_Ub2*SeqAggP;
            dwdp[233] = 1.0*Parkin_asyn_dam_Ub3*SeqAggP;
            dwdp[234] = 1.0*Parkin_asyn_dam_Ub4*SeqAggP;
            dwdp[235] = 1.0*Parkin_asyn_dam_Ub5*SeqAggP;
            dwdp[236] = 1.0*Parkin_asyn_dam_Ub6*SeqAggP;
            dwdp[237] = 1.0*Parkin_asyn_dam_Ub7*SeqAggP;
            dwdp[238] = 1.0*Parkin_asyn_dam_Ub8*SeqAggP;
            dwdp[239] = 1.0*Parkin_asyn_dam_Ub_DUB*SeqAggP;
            dwdp[240] = 1.0*Parkin_asyn_dam_Ub2_DUB*SeqAggP;
            dwdp[241] = 1.0*Parkin_asyn_dam_Ub3_DUB*SeqAggP;
            dwdp[242] = 1.0*Parkin_asyn_dam_Ub4_DUB*SeqAggP;
            dwdp[243] = 1.0*Parkin_asyn_dam_Ub5_DUB*SeqAggP;
            dwdp[244] = 1.0*Parkin_asyn_dam_Ub6_DUB*SeqAggP;
            dwdp[245] = 1.0*Parkin_asyn_dam_Ub7_DUB*SeqAggP;
            dwdp[246] = 1.0*Parkin_asyn_dam_Ub8_DUB*SeqAggP;
            dwdp[292] = 1.0*E3SUB_SUB_misfolded*SeqAggP;
            dwdp[293] = 1.0*E3SUB_SUB_misfolded_Ub*SeqAggP;
            dwdp[294] = 1.0*E3SUB_SUB_misfolded_Ub2*SeqAggP;
            dwdp[295] = 1.0*E3SUB_SUB_misfolded_Ub3*SeqAggP;
            dwdp[296] = 1.0*E3SUB_SUB_misfolded_Ub4*SeqAggP;
            dwdp[297] = 1.0*E3SUB_SUB_misfolded_Ub5*SeqAggP;
            dwdp[298] = 1.0*E3SUB_SUB_misfolded_Ub6*SeqAggP;
            dwdp[299] = 1.0*E3SUB_SUB_misfolded_Ub7*SeqAggP;
            dwdp[300] = 1.0*E3SUB_SUB_misfolded_Ub8*SeqAggP;
            dwdp[301] = 1.0*E3SUB_SUB_misfolded_Ub_UCHL1*SeqAggP;
            dwdp[302] = 1.0*E3SUB_SUB_misfolded_Ub2_UCHL1*SeqAggP;
            dwdp[303] = 1.0*E3SUB_SUB_misfolded_Ub3_UCHL1*SeqAggP;
            dwdp[304] = 1.0*E3SUB_SUB_misfolded_Ub4_UCHL1*SeqAggP;
            dwdp[305] = 1.0*E3SUB_SUB_misfolded_Ub5_UCHL1*SeqAggP;
            dwdp[306] = 1.0*E3SUB_SUB_misfolded_Ub6_UCHL1*SeqAggP;
            dwdp[307] = 1.0*E3SUB_SUB_misfolded_Ub7_UCHL1*SeqAggP;
            dwdp[308] = 1.0*E3SUB_SUB_misfolded_Ub8_UCHL1*SeqAggP;
            break;
        case 22:
            dwdp[78] = 1.0*AggP1*Proteasome;
            dwdp[79] = 1.0*AggP2*Proteasome;
            dwdp[80] = 1.0*AggP3*Proteasome;
            dwdp[81] = 1.0*AggP4*Proteasome;
            dwdp[82] = 1.0*AggP5*Proteasome;
            dwdp[201] = 1.0*AggA1*Proteasome;
            dwdp[202] = 1.0*AggA2*Proteasome;
            dwdp[203] = 1.0*AggA3*Proteasome;
            dwdp[204] = 1.0*AggA4*Proteasome;
            dwdp[205] = 1.0*AggA5*Proteasome;
            dwdp[223] = 1.0*AggD1*Proteasome;
            dwdp[224] = 1.0*AggD2*Proteasome;
            dwdp[225] = 1.0*AggD3*Proteasome;
            dwdp[226] = 1.0*AggD4*Proteasome;
            dwdp[227] = 1.0*AggD5*Proteasome;
            dwdp[262] = 1.0*AggU1*Proteasome;
            dwdp[263] = 1.0*AggU2*Proteasome;
            dwdp[264] = 1.0*AggU3*Proteasome;
            dwdp[265] = 1.0*AggU4*Proteasome;
            dwdp[266] = 1.0*AggU5*Proteasome;
            dwdp[285] = 1.0*AggS1*Proteasome;
            dwdp[286] = 1.0*AggS2*Proteasome;
            dwdp[287] = 1.0*AggS3*Proteasome;
            dwdp[288] = 1.0*AggS4*Proteasome;
            dwdp[289] = 1.0*AggS5*Proteasome;
            break;
        case 23:
            dwdp[18] = 1.0*DUB*E3_MisP_Ub;
            dwdp[19] = 1.0*DUB*E3_MisP_Ub2;
            dwdp[20] = 1.0*DUB*E3_MisP_Ub3;
            dwdp[21] = 1.0*DUB*E3_MisP_Ub4;
            dwdp[22] = 1.0*DUB*E3_MisP_Ub5;
            dwdp[23] = 1.0*DUB*E3_MisP_Ub6;
            dwdp[24] = 1.0*DUB*E3_MisP_Ub7;
            dwdp[25] = 1.0*DUB*E3_MisP_Ub8;
            break;
        case 24:
            dwdp[314] = 1.0*Source;
            break;
        case 25:
            dwdp[315] = 1.0*ROS;
            break;
        case 26:
            dwdp[0] = 1.0*Source;
            break;
        case 27:
            dwdp[1] = 1.0*Proteasome*Ub*kproteff;
            break;
        case 28:
            dwdp[2] = 1.0*pow(MisP, 6)/(1.0*pow(MisP, 6) + 11390625000000000000);
            break;
        case 29:
            dwdp[88] = 1.0*Source;
            break;
        case 30:
            dwdp[89] = 1.0*Proteasome*UCHL1;
            dwdp[93] = 1.0*Proteasome*UCHL1_damaged;
            break;
        case 31:
            dwdp[90] = 1.0*UCHL1_Proteasome*kproteff;
            dwdp[94] = 1.0*UCHL1_damaged_Proteasome*kproteff;
            break;
        case 32:
            dwdp[91] = 1.0*Lysosome*UCHL1;
            break;
        case 33:
            dwdp[92] = 1.0*ROS*UCHL1;
            break;
        case 34:
            dwdp[96] = 1.0*Lamp2a*UCHL1_damaged;
            break;
        case 35:
            dwdp[97] = 1.0*Lamp2a_UCHL1_damaged;
            break;
        case 36:
            dwdp[95] = 1.0*Lysosome*UCHL1_damaged;
            break;
        case 37:
            dwdp[98] = 1.0*UCHL1*Ub;
            break;
        case 38:
            dwdp[99] = 1.0*Ub_UCHL1;
            break;
        case 39:
            dwdp[121] = 1.0*E3SUB_SUB_misfolded_Ub8_UCHL1;
            dwdp[122] = 1.0*E3SUB_SUB_misfolded_Ub7_UCHL1;
            dwdp[123] = 1.0*E3SUB_SUB_misfolded_Ub6_UCHL1;
            dwdp[124] = 1.0*E3SUB_SUB_misfolded_Ub5_UCHL1;
            dwdp[125] = 1.0*E3SUB_SUB_misfolded_Ub4_UCHL1;
            dwdp[126] = 1.0*E3SUB_SUB_misfolded_Ub3_UCHL1;
            dwdp[127] = 1.0*E3SUB_SUB_misfolded_Ub2_UCHL1;
            dwdp[128] = 1.0*E3SUB_SUB_misfolded_Ub_UCHL1;
            break;
        case 40:
            dwdp[113] = 1.0*E3SUB_SUB_misfolded_Ub*UCHL1;
            dwdp[114] = 1.0*E3SUB_SUB_misfolded_Ub2*UCHL1;
            dwdp[115] = 1.0*E3SUB_SUB_misfolded_Ub3*UCHL1;
            dwdp[116] = 1.0*E3SUB_SUB_misfolded_Ub4*UCHL1;
            dwdp[117] = 1.0*E3SUB_SUB_misfolded_Ub5*UCHL1;
            dwdp[118] = 1.0*E3SUB_SUB_misfolded_Ub6*UCHL1;
            dwdp[119] = 1.0*E3SUB_SUB_misfolded_Ub7*UCHL1;
            dwdp[120] = 1.0*E3SUB_SUB_misfolded_Ub8*UCHL1;
            break;
        case 41:
            dwdp[100] = 1.0*Source;
            break;
        case 42:
            dwdp[101] = 1.0*ROS*SUB;
            break;
        case 43:
            dwdp[102] = 1.0*SUB_misfolded;
            break;
        case 44:
            dwdp[103] = 1.0*E3SUB*SUB_misfolded;
            break;
        case 45:
            dwdp[104] = 1.0*E3SUB_SUB_misfolded;
            break;
        case 46:
            dwdp[144] = 1.0*Source;
            break;
        case 47:
            dwdp[145] = 1.0*Proteasome*asyn;
            break;
        case 48:
            dwdp[146] = 1.0*asyn_Proteasome*kproteff;
            break;
        case 49:
            dwdp[147] = 1.0*Lamp2a*asyn;
            break;
        case 50:
            dwdp[148] = 1.0*asyn_Lamp2a;
            break;
        case 51:
            dwdp[149] = 1.0*ROS*asyn;
            break;
        case 52:
            dwdp[160] = 1.0*DUB*Parkin_asyn_dam_Ub8;
            dwdp[161] = 1.0*DUB*Parkin_asyn_dam_Ub7;
            dwdp[162] = 1.0*DUB*Parkin_asyn_dam_Ub6;
            dwdp[163] = 1.0*DUB*Parkin_asyn_dam_Ub5;
            dwdp[164] = 1.0*DUB*Parkin_asyn_dam_Ub4;
            dwdp[165] = 1.0*DUB*Parkin_asyn_dam_Ub3;
            dwdp[166] = 1.0*DUB*Parkin_asyn_dam_Ub2;
            dwdp[167] = 1.0*DUB*Parkin_asyn_dam_Ub;
            break;
        case 53:
            dwdp[150] = 1.0*Parkin*asyn_dam;
            break;
        case 54:
            dwdp[151] = 1.0*Parkin_asyn_dam;
            break;
        case 55:
            dwdp[191] = 0.5*asyn*(1.0*asyn - 1);
            break;
        case 56:
            dwdp[192] = 1.0*AggA1*asyn;
            dwdp[193] = 1.0*AggA2*asyn;
            dwdp[194] = 1.0*AggA3*asyn;
            dwdp[195] = 1.0*AggA4*asyn;
            dwdp[206] = 1.0*AggA5*asyn;
            break;
        case 57:
            dwdp[200] = 1.0*AggA1;
            break;
        case 58:
            dwdp[199] = 1.0*AggA2;
            break;
        case 59:
            dwdp[198] = 1.0*AggA3;
            break;
        case 60:
            dwdp[197] = 1.0*AggA4;
            break;
        case 61:
            dwdp[196] = 1.0*AggA5;
            break;
        case 62:
            dwdp[83] = 1.0*AggP1;
            dwdp[84] = 1.0*AggP2;
            dwdp[85] = 1.0*AggP3;
            dwdp[86] = 1.0*AggP4;
            dwdp[87] = 1.0*AggP5;
            dwdp[208] = 1.0*AggA1;
            dwdp[209] = 1.0*AggA2;
            dwdp[210] = 1.0*AggA3;
            dwdp[211] = 1.0*AggA4;
            dwdp[212] = 1.0*AggA5;
            dwdp[247] = 1.0*AggD1;
            dwdp[248] = 1.0*AggD2;
            dwdp[249] = 1.0*AggD3;
            dwdp[250] = 1.0*AggD4;
            dwdp[251] = 1.0*AggD5;
            dwdp[270] = 1.0*AggU1;
            dwdp[271] = 1.0*AggU2;
            dwdp[272] = 1.0*AggU3;
            dwdp[273] = 1.0*AggU4;
            dwdp[274] = 1.0*AggU5;
            dwdp[309] = 1.0*AggS1;
            dwdp[310] = 1.0*AggS2;
            dwdp[311] = 1.0*AggS3;
            dwdp[312] = 1.0*AggS4;
            dwdp[313] = 1.0*AggS5;
            break;
        case 63:
            dwdp[213] = 0.5*asyn_dam*(1.0*asyn_dam - 1);
            dwdp[252] = 0.5*UCHL1_damaged*(1.0*UCHL1_damaged - 1);
            break;
        case 64:
            dwdp[214] = 1.0*AggD1*asyn_dam;
            dwdp[215] = 1.0*AggD2*asyn_dam;
            dwdp[216] = 1.0*AggD3*asyn_dam;
            dwdp[217] = 1.0*AggD4*asyn_dam;
            dwdp[228] = 1.0*AggD5*asyn_dam;
            dwdp[253] = 1.0*AggU1*UCHL1_damaged;
            dwdp[254] = 1.0*AggU2*UCHL1_damaged;
            dwdp[255] = 1.0*AggU3*UCHL1_damaged;
            dwdp[256] = 1.0*AggU4*UCHL1_damaged;
            dwdp[267] = 1.0*AggU5*UCHL1_damaged;
            break;
        case 65:
            dwdp[222] = 1.0*AggD1;
            break;
        case 66:
            dwdp[221] = 1.0*AggD2;
            break;
        case 67:
            dwdp[220] = 1.0*AggD3;
            break;
        case 68:
            dwdp[219] = 1.0*AggD4;
            break;
        case 69:
            dwdp[218] = 1.0*AggD5;
            break;
        case 70:
            dwdp[261] = 1.0*AggU1;
            break;
        case 71:
            dwdp[260] = 1.0*AggU2;
            break;
        case 72:
            dwdp[259] = 1.0*AggU3;
            break;
        case 73:
            dwdp[258] = 1.0*AggU4;
            break;
        case 74:
            dwdp[257] = 1.0*AggU5;
            break;
        case 75:
            dwdp[275] = 0.5*SUB_misfolded*(1.0*SUB_misfolded - 1);
            break;
        case 76:
            dwdp[276] = 1.0*AggS1*SUB_misfolded;
            dwdp[277] = 1.0*AggS2*SUB_misfolded;
            dwdp[278] = 1.0*AggS3*SUB_misfolded;
            dwdp[279] = 1.0*AggS4*SUB_misfolded;
            dwdp[290] = 1.0*AggS5*SUB_misfolded;
            break;
        case 77:
            dwdp[284] = 1.0*AggS1;
            break;
        case 78:
            dwdp[283] = 1.0*AggS2;
            break;
        case 79:
            dwdp[282] = 1.0*AggS3;
            break;
        case 80:
            dwdp[281] = 1.0*AggS4;
            break;
        case 81:
            dwdp[280] = 1.0*AggS5;
            break;
        case 82:
            dwdp[1] = 1.0*Proteasome*Ub*kubd;
            dwdp[44] = 1.0*ATP*MisP_Ub4_Proteasome*kactProt/(1.0*ATP + 5000);
            dwdp[45] = 1.0*ATP*MisP_Ub5_Proteasome*kactProt/(1.0*ATP + 5000);
            dwdp[46] = 1.0*ATP*MisP_Ub6_Proteasome*kactProt/(1.0*ATP + 5000);
            dwdp[47] = 1.0*ATP*MisP_Ub7_Proteasome*kactProt/(1.0*ATP + 5000);
            dwdp[48] = 1.0*ATP*MisP_Ub8_Proteasome*kactProt/(1.0*ATP + 5000);
            dwdp[90] = 1.0*UCHL1_Proteasome*kdegProtUCHL1;
            dwdp[94] = 1.0*UCHL1_damaged_Proteasome*kdegProtUCHL1;
            dwdp[139] = 1.0*ATP*SUB_misfolded_Ub4_Proteasome*kactProt/(1.0*ATP + 5000);
            dwdp[140] = 1.0*ATP*SUB_misfolded_Ub5_Proteasome*kactProt/(1.0*ATP + 5000);
            dwdp[141] = 1.0*ATP*SUB_misfolded_Ub6_Proteasome*kactProt/(1.0*ATP + 5000);
            dwdp[142] = 1.0*ATP*SUB_misfolded_Ub7_Proteasome*kactProt/(1.0*ATP + 5000);
            dwdp[143] = 1.0*ATP*SUB_misfolded_Ub8_Proteasome*kactProt/(1.0*ATP + 5000);
            dwdp[146] = 1.0*asyn_Proteasome*kdegasynProt;
            dwdp[186] = 1.0*ATP*asyn_dam_Ub4_Proteasome*kactProt/(1.0*ATP + 5000);
            dwdp[187] = 1.0*ATP*asyn_dam_Ub5_Proteasome*kactProt/(1.0*ATP + 5000);
            dwdp[188] = 1.0*ATP*asyn_dam_Ub6_Proteasome*kactProt/(1.0*ATP + 5000);
            dwdp[189] = 1.0*ATP*asyn_dam_Ub7_Proteasome*kactProt/(1.0*ATP + 5000);
            dwdp[190] = 1.0*ATP*asyn_dam_Ub8_Proteasome*kactProt/(1.0*ATP + 5000);
            break;
    }
}