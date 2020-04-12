#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dxdotdw_model29_beuke30(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdotdw[0] = -76923.076923076922;
    dxdotdw[1] = -76923.076923076922;
    dxdotdw[2] = 76923.076923076922;
    dxdotdw[3] = 76923.076923076922;
    dxdotdw[4] = -76923.076923076922;
    dxdotdw[5] = -1250000.0;
    dxdotdw[6] = 1250000.0;
    dxdotdw[7] = 76923.076923076922;
    dxdotdw[8] = -76923.076923076922;
    dxdotdw[9] = -1250000.0;
    dxdotdw[10] = 1250000.0;
    dxdotdw[11] = -1250000.0;
    dxdotdw[12] = 76923.076923076922;
    dxdotdw[13] = -76923.076923076922;
    dxdotdw[14] = 76923.076923076922;
    dxdotdw[15] = -76923.076923076922;
    dxdotdw[16] = -1250000.0;
    dxdotdw[17] = -76923.076923076922;
    dxdotdw[18] = -1250000.0;
    dxdotdw[19] = 76923.076923076922;
    dxdotdw[20] = -1250000.0;
    dxdotdw[21] = -76923.076923076922;
    dxdotdw[22] = 1250000.0;
    dxdotdw[23] = -76923.076923076922;
    dxdotdw[24] = 76923.076923076922;
    dxdotdw[25] = -1250000.0;
    dxdotdw[26] = 1250000.0;
    dxdotdw[27] = 1250000.0;
    dxdotdw[28] = -1250000.0;
    dxdotdw[29] = 76923.076923076922;
    dxdotdw[30] = -76923.076923076922;
    dxdotdw[31] = 76923.076923076922;
    dxdotdw[32] = -76923.076923076922;
    dxdotdw[33] = 76923.076923076922;
    dxdotdw[34] = 76923.076923076922;
    dxdotdw[35] = -76923.076923076922;
    dxdotdw[36] = -386100.38610038609;
    dxdotdw[37] = -386100.38610038609;
    dxdotdw[38] = -1250000.0;
    dxdotdw[39] = 1250000.0;
    dxdotdw[40] = 1250000.0;
    dxdotdw[41] = -1250000.0;
    dxdotdw[42] = 386100.38610038609;
    dxdotdw[43] = -76923.076923076922;
    dxdotdw[44] = 2127659.5744680851;
    dxdotdw[45] = -2127659.5744680851;
    dxdotdw[46] = 4255319.1489361702;
    dxdotdw[47] = 386100.38610038609;
    dxdotdw[48] = 2857142.8571428573;
    dxdotdw[49] = -2857142.8571428573;
    dxdotdw[50] = 5714285.7142857146;
    dxdotdw[51] = 386100.38610038609;
    dxdotdw[52] = -386100.38610038609;
    dxdotdw[53] = 386100.38610038609;
    dxdotdw[54] = -386100.38610038609;
    dxdotdw[55] = -386100.38610038609;
    dxdotdw[56] = 76923.076923076922;
    dxdotdw[57] = -386100.38610038609;
    dxdotdw[58] = 76923.076923076922;
    dxdotdw[59] = -76923.076923076922;
    dxdotdw[60] = -2127659.5744680851;
    dxdotdw[61] = -2857142.8571428573;
    dxdotdw[62] = -76923.076923076922;
    dxdotdw[63] = 76923.076923076922;
    dxdotdw[64] = -76923.076923076922;
    dxdotdw[65] = 76923.076923076922;
    dxdotdw[66] = -76923.076923076922;
    dxdotdw[67] = 76923.076923076922;
    dxdotdw[68] = 76923.076923076922;
    dxdotdw[69] = -76923.076923076922;
    dxdotdw[70] = -1250000.0;
    dxdotdw[71] = -1250000.0;
    dxdotdw[72] = 1250000.0;
    dxdotdw[73] = -76923.076923076922;
    dxdotdw[74] = -1250000.0;
    dxdotdw[75] = 1250000.0;
    dxdotdw[76] = -1250000.0;
    dxdotdw[77] = 76923.076923076922;
    dxdotdw[78] = -76923.076923076922;
    dxdotdw[79] = 1250000.0;
    dxdotdw[80] = -1250000.0;
    dxdotdw[81] = -76923.076923076922;
    dxdotdw[82] = 1250000.0;
    dxdotdw[83] = 1250000.0;
    dxdotdw[84] = -1250000.0;
    dxdotdw[85] = 76923.076923076922;
    dxdotdw[86] = 1250000.0;
    dxdotdw[87] = -76923.076923076922;
    dxdotdw[88] = -76923.076923076922;
    dxdotdw[89] = 76923.076923076922;
    dxdotdw[90] = 76923.076923076922;
    dxdotdw[91] = 76923.076923076922;
    dxdotdw[92] = -76923.076923076922;
    dxdotdw[93] = 1250000.0;
    dxdotdw[94] = -1250000.0;
    dxdotdw[95] = -1250000.0;
    dxdotdw[96] = -76923.076923076922;
    dxdotdw[97] = 1250000.0;
}