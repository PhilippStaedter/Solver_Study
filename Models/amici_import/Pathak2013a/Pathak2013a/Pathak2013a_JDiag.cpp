#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_Pathak2013a(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2;
    JDiag[1] = -1.0*dwdx3 - 1.0*dwdx4 - 1.0*dwdx5;
    JDiag[2] = 1.0*dwdx6 - 1.0*dwdx7;
    JDiag[3] = -1.0*dwdx10 + 1.0*dwdx8 + 1.0*dwdx9;
    JDiag[4] = 1.0*dwdx11 + 1.0*dwdx12 - 1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx15;
    JDiag[5] = 1.0*dwdx16 - 1.0*dwdx17 - 1.0*dwdx18;
    JDiag[6] = -1.0*dwdx19 + 1.0*dwdx20 + 1.0*dwdx21 + 1.0*dwdx22 + 1.0*dwdx23;
    JDiag[7] = 1.0*dwdx24 - 1.0*dwdx25 - 1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28 - 1.0*dwdx29 - 1.0*dwdx30;
    JDiag[8] = 1.0*dwdx31 + 1.0*dwdx32 - 1.0*dwdx33 - 1.0*dwdx34 - 1.0*dwdx35;
    JDiag[9] = 1.0*dwdx36;
    JDiag[10] = 1.0*dwdx37 + 1.0*dwdx38 - 1.0*dwdx39;
    JDiag[11] = 1.0*dwdx40;
    JDiag[12] = 1.0*dwdx41 + 1.0*dwdx42;
    JDiag[13] = -1.0*dwdx43 + 1.0*dwdx44;
    JDiag[14] = 1.0*dwdx45 - 1.0*dwdx46 - 1.0*dwdx47 - 1.0*dwdx48 - 1.0*dwdx49 - 1.0*dwdx50;
    JDiag[15] = 1.0*dwdx51 + 1.0*dwdx52 - 1.0*dwdx53 - 1.0*dwdx54;
    JDiag[16] = 1.0*dwdx55 + 1.0*dwdx56 - 1.0*dwdx57 - 1.0*dwdx58 - 1.0*dwdx59;
    JDiag[17] = 1.0*dwdx60 + 1.0*dwdx61 - 1.0*dwdx62 - 1.0*dwdx63;
    JDiag[18] = 1.0*dwdx64 + 1.0*dwdx65;
    JDiag[19] = -1.0*dwdx66 + 1.0*dwdx67;
    JDiag[20] = 1.0*dwdx68 - 1.0*dwdx69 - 1.0*dwdx70 - 1.0*dwdx71 - 1.0*dwdx72 - 1.0*dwdx73 - 1.0*dwdx74 - 1.0*dwdx75 - 1.0*dwdx76 - 1.0*dwdx77 - 1.0*dwdx78 - 1.0*dwdx79;
    JDiag[21] = 1.0*dwdx80 + 1.0*dwdx81 + 1.0*dwdx82 - 1.0*dwdx83;
    JDiag[22] = 1.0*dwdx84 + 1.0*dwdx85 + 1.0*dwdx86 - 1.0*dwdx87 - 1.0*dwdx88 - 1.0*dwdx89 - 1.0*dwdx90 - 1.0*dwdx91;
    JDiag[23] = 1.0*dwdx92 + 1.0*dwdx93 - 1.0*dwdx94 - 1.0*dwdx95 - 1.0*dwdx96 - 1.0*dwdx97;
    JDiag[24] = 1.0*dwdx100 - 1.0*dwdx101 - 1.0*dwdx102 - 1.0*dwdx103 - 1.0*dwdx104 - 1.0*dwdx105 + 1.0*dwdx98 + 1.0*dwdx99;
    JDiag[25] = -1.0*dwdx106 + 1.0*dwdx107 + 1.0*dwdx108;
    JDiag[26] = 1.0*dwdx109 - 1.0*dwdx110;
    JDiag[27] = -1.0*dwdx111 + 1.0*dwdx112;
    JDiag[28] = 1.0*dwdx113 - 1.0*dwdx114;
    JDiag[29] = -1.0*dwdx115 + 1.0*dwdx116 + 1.0*dwdx117 + 1.0*dwdx118;
    JDiag[30] = 1.0*dwdx119 - 1.0*dwdx120;
    JDiag[31] = -1.0*dwdx121 + 1.0*dwdx122;
    JDiag[32] = 1.0*dwdx123 - 1.0*dwdx124;
    JDiag[33] = -1.0*dwdx125 + 1.0*dwdx126;
    JDiag[34] = 1.0*dwdx127 - 1.0*dwdx128;
    JDiag[35] = -1.0*dwdx129 + 1.0*dwdx130;
    JDiag[36] = 1.0*dwdx131 - 1.0*dwdx132;
    JDiag[37] = -1.0*dwdx133 + 1.0*dwdx134;
    JDiag[38] = 1.0*dwdx135 - 1.0*dwdx136;
    JDiag[39] = -1.0*dwdx137 + 1.0*dwdx138 + 1.0*dwdx139;
    JDiag[40] = 1.0*dwdx140 - 1.0*dwdx141;
    JDiag[41] = -1.0*dwdx142 + 1.0*dwdx143 + 1.0*dwdx144;
    JDiag[42] = 1.0*dwdx145 - 1.0*dwdx146;
    JDiag[43] = -1.0*dwdx147 + 1.0*dwdx148 + 1.0*dwdx149;
    JDiag[44] = 1.0*dwdx150 - 1.0*dwdx151;
    JDiag[45] = 1.0*dwdx152 - 1.0*dwdx153;
    JDiag[46] = 1.0*dwdx154 - 1.0*dwdx155;
    JDiag[47] = 1.0*dwdx156 - 1.0*dwdx157;
    JDiag[48] = 1.0*dwdx158 - 1.0*dwdx159;
    JDiag[49] = 1.0*dwdx160 + 1.0*dwdx161 + 1.0*dwdx162 + 1.0*dwdx163 + 1.0*dwdx164 + 1.0*dwdx165 + 1.0*dwdx166 + 1.0*dwdx167 + 1.0*dwdx168 + 1.0*dwdx169 + 1.0*dwdx170 + 1.0*dwdx171 + 1.0*dwdx172 + 1.0*dwdx173;
    JDiag[50] = 1.0*dwdx174;
    JDiag[51] = 1.0*dwdx175;
}