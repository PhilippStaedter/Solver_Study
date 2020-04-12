#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void sigmay_Lai2014(realtype *sigmay, const realtype t, const realtype *p, const realtype *k){
    sigmay[0] = 1.0;
    sigmay[1] = 1.0;
    sigmay[2] = 1.0;
    sigmay[3] = 1.0;
    sigmay[4] = 1.0;
    sigmay[5] = 1.0;
    sigmay[6] = 1.0;
    sigmay[7] = 1.0;
    sigmay[8] = 1.0;
    sigmay[9] = 1.0;
    sigmay[10] = 1.0;
    sigmay[11] = 1.0;
    sigmay[12] = 1.0;
    sigmay[13] = 1.0;
    sigmay[14] = 1.0;
    sigmay[15] = 1.0;
    sigmay[16] = 1.0;
    sigmay[17] = 1.0;
    sigmay[18] = 1.0;
    sigmay[19] = 1.0;
    sigmay[20] = 1.0;
    sigmay[21] = 1.0;
    sigmay[22] = 1.0;
    sigmay[23] = 1.0;
    sigmay[24] = 1.0;
    sigmay[25] = 1.0;
    sigmay[26] = 1.0;
    sigmay[27] = 1.0;
    sigmay[28] = 1.0;
    sigmay[29] = 1.0;
    sigmay[30] = 1.0;
    sigmay[31] = 1.0;
    sigmay[32] = 1.0;
    sigmay[33] = 1.0;
    sigmay[34] = 1.0;
    sigmay[35] = 1.0;
    sigmay[36] = 1.0;
    sigmay[37] = 1.0;
    sigmay[38] = 1.0;
    sigmay[39] = 1.0;
    sigmay[40] = 1.0;
    sigmay[41] = 1.0;
    sigmay[42] = 1.0;
    sigmay[43] = 1.0;
    sigmay[44] = 1.0;
    sigmay[45] = 1.0;
    sigmay[46] = 1.0;
    sigmay[47] = 1.0;
    sigmay[48] = 1.0;
    sigmay[49] = 1.0;
    sigmay[50] = 1.0;
    sigmay[51] = 1.0;
    sigmay[52] = 1.0;
    sigmay[53] = 1.0;
    sigmay[54] = 1.0;
    sigmay[55] = 1.0;
    sigmay[56] = 1.0;
    sigmay[57] = 1.0;
    sigmay[58] = 1.0;
    sigmay[59] = 1.0;
    sigmay[60] = 1.0;
    sigmay[61] = 1.0;
    sigmay[62] = 1.0;
    sigmay[63] = 1.0;
    sigmay[64] = 1.0;
    sigmay[65] = 1.0;
    sigmay[66] = 1.0;
    sigmay[67] = 1.0;
    sigmay[68] = 1.0;
    sigmay[69] = 1.0;
    sigmay[70] = 1.0;
    sigmay[71] = 1.0;
    sigmay[72] = 1.0;
    sigmay[73] = 1.0;
    sigmay[74] = 1.0;
    sigmay[75] = 1.0;
    sigmay[76] = 1.0;
    sigmay[77] = 1.0;
    sigmay[78] = 1.0;
    sigmay[79] = 1.0;
    sigmay[80] = 1.0;
    sigmay[81] = 1.0;
    sigmay[82] = 1.0;
    sigmay[83] = 1.0;
    sigmay[84] = 1.0;
    sigmay[85] = 1.0;
    sigmay[86] = 1.0;
    sigmay[87] = 1.0;
    sigmay[88] = 1.0;
    sigmay[89] = 1.0;
    sigmay[90] = 1.0;
    sigmay[91] = 1.0;
    sigmay[92] = 1.0;
    sigmay[93] = 1.0;
    sigmay[94] = 1.0;
    sigmay[95] = 1.0;
    sigmay[96] = 1.0;
    sigmay[97] = 1.0;
    sigmay[98] = 1.0;
    sigmay[99] = 1.0;
    sigmay[100] = 1.0;
    sigmay[101] = 1.0;
    sigmay[102] = 1.0;
    sigmay[103] = 1.0;
    sigmay[104] = 1.0;
    sigmay[105] = 1.0;
    sigmay[106] = 1.0;
    sigmay[107] = 1.0;
    sigmay[108] = 1.0;
    sigmay[109] = 1.0;
    sigmay[110] = 1.0;
    sigmay[111] = 1.0;
    sigmay[112] = 1.0;
    sigmay[113] = 1.0;
    sigmay[114] = 1.0;
    sigmay[115] = 1.0;
    sigmay[116] = 1.0;
    sigmay[117] = 1.0;
    sigmay[118] = 1.0;
    sigmay[119] = 1.0;
    sigmay[120] = 1.0;
    sigmay[121] = 1.0;
    sigmay[122] = 1.0;
    sigmay[123] = 1.0;
    sigmay[124] = 1.0;
    sigmay[125] = 1.0;
    sigmay[126] = 1.0;
    sigmay[127] = 1.0;
    sigmay[128] = 1.0;
    sigmay[129] = 1.0;
    sigmay[130] = 1.0;
    sigmay[131] = 1.0;
    sigmay[132] = 1.0;
    sigmay[133] = 1.0;
    sigmay[134] = 1.0;
    sigmay[135] = 1.0;
    sigmay[136] = 1.0;
    sigmay[137] = 1.0;
    sigmay[138] = 1.0;
    sigmay[139] = 1.0;
    sigmay[140] = 1.0;
    sigmay[141] = 1.0;
    sigmay[142] = 1.0;
    sigmay[143] = 1.0;
    sigmay[144] = 1.0;
    sigmay[145] = 1.0;
    sigmay[146] = 1.0;
    sigmay[147] = 1.0;
    sigmay[148] = 1.0;
    sigmay[149] = 1.0;
    sigmay[150] = 1.0;
    sigmay[151] = 1.0;
    sigmay[152] = 1.0;
    sigmay[153] = 1.0;
    sigmay[154] = 1.0;
    sigmay[155] = 1.0;
    sigmay[156] = 1.0;
    sigmay[157] = 1.0;
    sigmay[158] = 1.0;
    sigmay[159] = 1.0;
    sigmay[160] = 1.0;
    sigmay[161] = 1.0;
    sigmay[162] = 1.0;
    sigmay[163] = 1.0;
    sigmay[164] = 1.0;
    sigmay[165] = 1.0;
    sigmay[166] = 1.0;
    sigmay[167] = 1.0;
    sigmay[168] = 1.0;
    sigmay[169] = 1.0;
    sigmay[170] = 1.0;
    sigmay[171] = 1.0;
    sigmay[172] = 1.0;
    sigmay[173] = 1.0;
    sigmay[174] = 1.0;
    sigmay[175] = 1.0;
    sigmay[176] = 1.0;
    sigmay[177] = 1.0;
    sigmay[178] = 1.0;
    sigmay[179] = 1.0;
    sigmay[180] = 1.0;
    sigmay[181] = 1.0;
    sigmay[182] = 1.0;
    sigmay[183] = 1.0;
    sigmay[184] = 1.0;
    sigmay[185] = 1.0;
    sigmay[186] = 1.0;
    sigmay[187] = 1.0;
    sigmay[188] = 1.0;
    sigmay[189] = 1.0;
    sigmay[190] = 1.0;
    sigmay[191] = 1.0;
    sigmay[192] = 1.0;
    sigmay[193] = 1.0;
    sigmay[194] = 1.0;
}