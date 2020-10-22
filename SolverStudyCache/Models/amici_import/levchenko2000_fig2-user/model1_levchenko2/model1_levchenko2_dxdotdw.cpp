#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dxdotdw_model1_levchenko2(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdotdw[0] = -1.0;
    dxdotdw[1] = -1.0;
    dxdotdw[2] = 1.0;
    dxdotdw[3] = -1.0;
    dxdotdw[4] = -1.0;
    dxdotdw[5] = 1.0;
    dxdotdw[6] = 1.0;
    dxdotdw[7] = 1.0;
    dxdotdw[8] = -1.0;
    dxdotdw[9] = 1.0;
    dxdotdw[10] = 1.0;
    dxdotdw[11] = -1.0;
    dxdotdw[12] = -1.0;
    dxdotdw[13] = 1.0;
    dxdotdw[14] = -1.0;
    dxdotdw[15] = 1.0;
    dxdotdw[16] = -1.0;
    dxdotdw[17] = 1.0;
    dxdotdw[18] = -1.0;
    dxdotdw[19] = 1.0;
    dxdotdw[20] = 1.0;
    dxdotdw[21] = -1.0;
    dxdotdw[22] = -1.0;
    dxdotdw[23] = 1.0;
    dxdotdw[24] = 1.0;
    dxdotdw[25] = 1.0;
    dxdotdw[26] = -1.0;
    dxdotdw[27] = 1.0;
    dxdotdw[28] = 1.0;
    dxdotdw[29] = -1.0;
    dxdotdw[30] = -1.0;
    dxdotdw[31] = 1.0;
    dxdotdw[32] = -1.0;
    dxdotdw[33] = 1.0;
    dxdotdw[34] = 1.0;
    dxdotdw[35] = -1.0;
    dxdotdw[36] = 1.0;
    dxdotdw[37] = -1.0;
    dxdotdw[38] = 1.0;
    dxdotdw[39] = -1.0;
    dxdotdw[40] = 1.0;
    dxdotdw[41] = 1.0;
    dxdotdw[42] = -1.0;
    dxdotdw[43] = -1.0;
    dxdotdw[44] = 1.0;
    dxdotdw[45] = 1.0;
    dxdotdw[46] = 1.0;
    dxdotdw[47] = -1.0;
    dxdotdw[48] = 1.0;
    dxdotdw[49] = 1.0;
    dxdotdw[50] = -1.0;
    dxdotdw[51] = -1.0;
    dxdotdw[52] = 1.0;
    dxdotdw[53] = -1.0;
    dxdotdw[54] = 1.0;
    dxdotdw[55] = -1.0;
    dxdotdw[56] = 1.0;
    dxdotdw[57] = -1.0;
    dxdotdw[58] = 1.0;
    dxdotdw[59] = 1.0;
    dxdotdw[60] = -1.0;
    dxdotdw[61] = -1.0;
    dxdotdw[62] = 1.0;
    dxdotdw[63] = 1.0;
    dxdotdw[64] = 1.0;
    dxdotdw[65] = -1.0;
    dxdotdw[66] = 1.0;
    dxdotdw[67] = -1.0;
    dxdotdw[68] = 1.0;
    dxdotdw[69] = 1.0;
    dxdotdw[70] = 1.0;
    dxdotdw[71] = -1.0;
    dxdotdw[72] = -1.0;
    dxdotdw[73] = -1.0;
    dxdotdw[74] = 1.0;
    dxdotdw[75] = 1.0;
    dxdotdw[76] = 1.0;
    dxdotdw[77] = -1.0;
    dxdotdw[78] = 1.0;
    dxdotdw[79] = 1.0;
    dxdotdw[80] = -1.0;
    dxdotdw[81] = -1.0;
    dxdotdw[82] = 1.0;
    dxdotdw[83] = -1.0;
    dxdotdw[84] = 1.0;
    dxdotdw[85] = -1.0;
    dxdotdw[86] = 1.0;
    dxdotdw[87] = -1.0;
    dxdotdw[88] = 1.0;
    dxdotdw[89] = 1.0;
    dxdotdw[90] = -1.0;
    dxdotdw[91] = 1.0;
    dxdotdw[92] = -1.0;
    dxdotdw[93] = 1.0;
    dxdotdw[94] = -1.0;
    dxdotdw[95] = 1.0;
    dxdotdw[96] = -1.0;
    dxdotdw[97] = 1.0;
    dxdotdw[98] = -1.0;
    dxdotdw[99] = 1.0;
    dxdotdw[100] = -1.0;
    dxdotdw[101] = 1.0;
    dxdotdw[102] = -1.0;
    dxdotdw[103] = 1.0;
    dxdotdw[104] = -1.0;
    dxdotdw[105] = 1.0;
    dxdotdw[106] = -1.0;
    dxdotdw[107] = 1.0;
    dxdotdw[108] = 1.0;
    dxdotdw[109] = -1.0;
    dxdotdw[110] = 1.0;
    dxdotdw[111] = 1.0;
    dxdotdw[112] = -1.0;
    dxdotdw[113] = 1.0;
}