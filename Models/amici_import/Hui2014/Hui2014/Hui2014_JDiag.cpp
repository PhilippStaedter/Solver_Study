#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_Hui2014(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0;
    JDiag[1] = -1.0*dwdx2;
    JDiag[3] = -1.0*dwdx5 - 1.0*dwdx6;
    JDiag[4] = -1.0*dwdx7 - 1.0*dwdx8;
    JDiag[5] = -1.0*dwdx10 - 1.0*dwdx11 - 2.0*dwdx9;
    JDiag[6] = -1.0*dwdx12 - 1.0*dwdx13;
    JDiag[7] = -1.0*dwdx15 - 1.0*dwdx16 - 1.0*dwdx17;
    JDiag[8] = -1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20;
    JDiag[9] = -1.0*dwdx21 - 1.0*dwdx22;
    JDiag[10] = -1.0*dwdx23 - 1.0*dwdx24;
    JDiag[11] = -1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28 - 1.0*dwdx29 - 1.0*dwdx30 - 1.0*dwdx31;
    JDiag[12] = -1.0*dwdx33 - 1.0*dwdx34;
    JDiag[13] = -1.0*dwdx35 - 1.0*dwdx36;
    JDiag[14] = -1.0*dwdx38 - 1.0*dwdx39 - 1.0*dwdx40 - 1.0*dwdx41;
    JDiag[15] = -1.0*dwdx43 - 1.0*dwdx44;
    JDiag[16] = -1.0*dwdx45 - 1.0*dwdx46 - 1.0*dwdx47;
    JDiag[17] = -1.0*dwdx50 - 1.0*dwdx51 - 1.0*dwdx52;
    JDiag[18] = -1.0*dwdx53;
    JDiag[19] = -1.0*dwdx55;
    JDiag[20] = -1.0*dwdx57;
    JDiag[21] = -1.0*dwdx58 - 1.0*dwdx59;
    JDiag[22] = -1.0*dwdx61;
    JDiag[23] = -1.0*dwdx66 - 1.0*dwdx68;
    JDiag[24] = -1.0*dwdx69;
    JDiag[25] = -1.0*dwdx70;
    JDiag[26] = -1.0*dwdx72;
    JDiag[28] = -1.0*dwdx75 - 1.0*dwdx76;
    JDiag[29] = -1.0*dwdx81;
    JDiag[30] = -1.0*dwdx82 - 1.0*dwdx83;
    JDiag[31] = -1.0*dwdx85;
    JDiag[32] = -1.0*dwdx88;
    JDiag[33] = -1.0*dwdx89;
    JDiag[34] = -1.0*dwdx91;
    JDiag[35] = -1.0*dwdx93 - 1.0*dwdx96;
    JDiag[36] = -1.0*dwdx99;
    JDiag[37] = -1.0*dwdx101;
    JDiag[38] = -1.0*dwdx102;
    JDiag[39] = -1.0*dwdx103 - 1.0*dwdx104 - 1.0*dwdx105;
    JDiag[40] = -1.0*dwdx106;
    JDiag[41] = -1.0*dwdx108;
    JDiag[42] = -1.0*dwdx109 - 1.0*dwdx110;
    JDiag[43] = -1.0*dwdx111;
    JDiag[44] = -1.0*dwdx116 - 1.0*dwdx117;
    JDiag[45] = -1.0*dwdx118 - 1.0*dwdx120 - 1.0*dwdx121;
    JDiag[46] = -1.0*dwdx122;
    JDiag[47] = -1.0*dwdx124 - 1.0*dwdx125;
    JDiag[48] = -1.0*dwdx126;
    JDiag[49] = -1.0*dwdx130;
    JDiag[50] = -1.0*dwdx132 - 1.0*dwdx133 - 1.0*dwdx134;
    JDiag[51] = -1.0*dwdx135 - 1.0*dwdx137;
    JDiag[52] = -1.0*dwdx138 - 1.0*dwdx139;
    JDiag[53] = -1.0*dwdx140 - 1.0*dwdx141;
    JDiag[54] = -1.0*dwdx143 - 1.0*dwdx144;
    JDiag[56] = -1.0*dwdx145;
    JDiag[57] = -1.0*dwdx146;
    JDiag[59] = -1.0*dwdx147 - 1.0*dwdx148;
    JDiag[60] = -1.0*dwdx149;
    JDiag[61] = -1.0*dwdx151 - 1.0*dwdx152;
}