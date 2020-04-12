#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Qi2013a(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = NatP;
    y[1] = MisP;
    y[2] = Ub;
    y[3] = E1;
    y[4] = E2;
    y[5] = E3;
    y[6] = E3_MisP;
    y[7] = DUB;
    y[8] = Proteasome;
    y[9] = ROS;
    y[10] = E1_Ub;
    y[11] = E2_Ub;
    y[12] = E3_MisP_Ub;
    y[13] = E3_MisP_Ub2;
    y[14] = E3_MisP_Ub3;
    y[15] = E3_MisP_Ub4;
    y[16] = E3_MisP_Ub5;
    y[17] = E3_MisP_Ub6;
    y[18] = E3_MisP_Ub7;
    y[19] = E3_MisP_Ub8;
    y[20] = MisP_Ub4_Proteasome;
    y[21] = MisP_Ub5_Proteasome;
    y[22] = MisP_Ub6_Proteasome;
    y[23] = MisP_Ub7_Proteasome;
    y[24] = MisP_Ub8_Proteasome;
    y[25] = E3_MisP_Ub_DUB;
    y[26] = E3_MisP_Ub2_DUB;
    y[27] = E3_MisP_Ub3_DUB;
    y[28] = E3_MisP_Ub4_DUB;
    y[29] = E3_MisP_Ub5_DUB;
    y[30] = E3_MisP_Ub6_DUB;
    y[31] = E3_MisP_Ub7_DUB;
    y[32] = E3_MisP_Ub8_DUB;
    y[33] = AggP1;
    y[34] = AggP2;
    y[35] = AggP3;
    y[36] = AggP4;
    y[37] = AggP5;
    y[38] = SeqAggP;
    y[39] = AggP_Proteasome;
    y[40] = ATP;
    y[41] = ADP;
    y[42] = AMP;
    y[43] = UCHL1;
    y[44] = UCHL1_Proteasome;
    y[45] = UCHL1_damaged_Proteasome;
    y[46] = Lysosome;
    y[47] = UCHL1_damaged;
    y[48] = Lamp2a;
    y[49] = Lamp2a_UCHL1_damaged;
    y[50] = Ub_UCHL1;
    y[51] = SUB;
    y[52] = SUB_misfolded;
    y[53] = E3SUB;
    y[54] = E3SUB_SUB_misfolded;
    y[55] = E3SUB_SUB_misfolded_Ub;
    y[56] = E3SUB_SUB_misfolded_Ub2;
    y[57] = E3SUB_SUB_misfolded_Ub3;
    y[58] = E3SUB_SUB_misfolded_Ub4;
    y[59] = E3SUB_SUB_misfolded_Ub5;
    y[60] = E3SUB_SUB_misfolded_Ub6;
    y[61] = E3SUB_SUB_misfolded_Ub7;
    y[62] = E3SUB_SUB_misfolded_Ub8;
    y[63] = E3SUB_SUB_misfolded_Ub_UCHL1;
    y[64] = E3SUB_SUB_misfolded_Ub2_UCHL1;
    y[65] = E3SUB_SUB_misfolded_Ub3_UCHL1;
    y[66] = E3SUB_SUB_misfolded_Ub4_UCHL1;
    y[67] = E3SUB_SUB_misfolded_Ub5_UCHL1;
    y[68] = E3SUB_SUB_misfolded_Ub6_UCHL1;
    y[69] = E3SUB_SUB_misfolded_Ub7_UCHL1;
    y[70] = E3SUB_SUB_misfolded_Ub8_UCHL1;
    y[71] = SUB_misfolded_Ub4_Proteasome;
    y[72] = SUB_misfolded_Ub5_Proteasome;
    y[73] = SUB_misfolded_Ub6_Proteasome;
    y[74] = SUB_misfolded_Ub7_Proteasome;
    y[75] = SUB_misfolded_Ub8_Proteasome;
    y[76] = asyn;
    y[77] = asyn_Proteasome;
    y[78] = asyn_Lamp2a;
    y[79] = asyn_dam;
    y[80] = Parkin;
    y[81] = Parkin_asyn_dam;
    y[82] = Parkin_asyn_dam_Ub;
    y[83] = Parkin_asyn_dam_Ub2;
    y[84] = Parkin_asyn_dam_Ub3;
    y[85] = Parkin_asyn_dam_Ub4;
    y[86] = Parkin_asyn_dam_Ub5;
    y[87] = Parkin_asyn_dam_Ub6;
    y[88] = Parkin_asyn_dam_Ub7;
    y[89] = Parkin_asyn_dam_Ub8;
    y[90] = Parkin_asyn_dam_Ub_DUB;
    y[91] = Parkin_asyn_dam_Ub2_DUB;
    y[92] = Parkin_asyn_dam_Ub3_DUB;
    y[93] = Parkin_asyn_dam_Ub4_DUB;
    y[94] = Parkin_asyn_dam_Ub5_DUB;
    y[95] = Parkin_asyn_dam_Ub6_DUB;
    y[96] = Parkin_asyn_dam_Ub7_DUB;
    y[97] = Parkin_asyn_dam_Ub8_DUB;
    y[98] = asyn_dam_Ub4_Proteasome;
    y[99] = asyn_dam_Ub5_Proteasome;
    y[100] = asyn_dam_Ub6_Proteasome;
    y[101] = asyn_dam_Ub7_Proteasome;
    y[102] = asyn_dam_Ub8_Proteasome;
    y[103] = AggA1;
    y[104] = AggA2;
    y[105] = AggA3;
    y[106] = AggA4;
    y[107] = AggA5;
    y[108] = AggD1;
    y[109] = AggD2;
    y[110] = AggD3;
    y[111] = AggD4;
    y[112] = AggD5;
    y[113] = AggU1;
    y[114] = AggU2;
    y[115] = AggU3;
    y[116] = AggU4;
    y[117] = AggU5;
    y[118] = AggS1;
    y[119] = AggS2;
    y[120] = AggS3;
    y[121] = AggS4;
    y[122] = AggS5;
    y[123] = Source;
    y[124] = Sink;
    y[125] = aggasyn;
    y[126] = aggasyndam;
    y[127] = aggParkin;
    y[128] = aggUb;
    y[129] = aggE3;
    y[130] = aggDUB;
    y[131] = aggMisP;
    y[132] = aggUchl1;
    y[133] = aggUchl1dam;
    y[134] = aggSUB;
    y[135] = upregUb;
}