#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model3_levchenko2(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*MAPK*on2 + 1.0*MEK*on1;
    dwdx[1] = 1.0*MEK*on1;
    dwdx[2] = 1.0*MAPK*on2;
    dwdx[3] = 1.0*MAPK*on2;
    dwdx[4] = 1.0*MEK*on1;
    dwdx[5] = of1;
    dwdx[6] = 1.0*MAPK*on2 + 1.0*RAFp*kr1 + of1;
    dwdx[7] = 1.0*RAFp*kr1;
    dwdx[8] = 1.0*MAPK*on2;
    dwdx[9] = 1.0*MAPK*on2;
    dwdx[10] = of1;
    dwdx[11] = of3;
    dwdx[12] = 1.0*MAPK*on2 + of3;
    dwdx[13] = 1.0*MAPK*on2;
    dwdx[14] = 1.0*MAPK*on2;
    dwdx[15] = of3;
    dwdx[16] = of2;
    dwdx[17] = 1.0*MEK*on1 + of2;
    dwdx[18] = 1.0*MEK*on1;
    dwdx[19] = of2;
    dwdx[20] = 1.0*MEK*on1;
    dwdx[21] = of4;
    dwdx[22] = 1.0*MEK*on1 + of4;
    dwdx[23] = 1.0*MEK*on1;
    dwdx[24] = of4;
    dwdx[25] = 1.0*MEK*on1;
    dwdx[26] = of2;
    dwdx[27] = of1;
    dwdx[28] = 1.0*RAFp*kr1 + of1 + of2;
    dwdx[29] = 1.0*RAFp*kr1;
    dwdx[30] = of2;
    dwdx[31] = of1;
    dwdx[32] = of2;
    dwdx[33] = of3;
    dwdx[34] = kr2 + of2 + of3;
    dwdx[35] = kr2;
    dwdx[36] = of2;
    dwdx[37] = of3;
    dwdx[38] = of4;
    dwdx[39] = of3;
    dwdx[40] = of3 + of4;
    dwdx[41] = of4;
    dwdx[42] = of3;
    dwdx[43] = of4;
    dwdx[44] = of1;
    dwdx[45] = 1.0*RAFp*kr1;
    dwdx[46] = 1.0*RAFp*kr1 + of1 + of4;
    dwdx[47] = of4;
    dwdx[48] = of1;
    dwdx[49] = 1.0*MEKpp*Reaction19_a7;
    dwdx[50] = 1.0*C1*on2;
    dwdx[51] = 1.0*C2*on2;
    dwdx[52] = 1.0*C3*on2;
    dwdx[53] = 1.0*C1*on2;
    dwdx[54] = 1.0*C2*on2;
    dwdx[55] = 1.0*C3*on2;
    dwdx[56] = 1.0*on2*(C1 + C2 + C3);
    dwdx[57] = 1.0*Reaction20_d7;
    dwdx[58] = 1.0*Reaction21_k7;
    dwdx[59] = 1.0*MAPKp*Reaction22_a8;
    dwdx[60] = 1.0*MAPKpp*Reaction28_a10;
    dwdx[61] = 1.0*MAPKPH*Reaction22_a8;
    dwdx[62] = 1.0*MEKpp*Reaction25_a9;
    dwdx[63] = 1.0*Reaction23_d8;
    dwdx[64] = 1.0*Reaction24_k8;
    dwdx[65] = 1.0*Reaction26_d9;
    dwdx[66] = 1.0*Reaction27_k9;
    dwdx[67] = 1.0*MAPKPH*Reaction28_a10;
    dwdx[68] = 1.0*Reaction29_d10;
    dwdx[69] = 1.0*Reaction30_k10;
    dwdx[70] = 1.0*RAFp*Reaction7_a3;
    dwdx[71] = 1.0*C1*on1;
    dwdx[72] = 1.0*C1*on1;
    dwdx[73] = 1.0*C4*on1;
    dwdx[74] = 1.0*C5*on1;
    dwdx[75] = 1.0*C4*on1;
    dwdx[76] = 1.0*C5*on1;
    dwdx[77] = 1.0*on1*(C1 + C4 + C5);
    dwdx[78] = 1.0*MEKp*Reaction10_a4;
    dwdx[79] = 1.0*MEKpp*Reaction16_a6;
    dwdx[80] = 1.0*Reaction8_d3;
    dwdx[81] = 1.0*Reaction9_k3;
    dwdx[82] = 1.0*MEKPH*Reaction10_a4;
    dwdx[83] = 1.0*RAFp*Reaction13_a5;
    dwdx[84] = 1.0*Reaction11_d4;
    dwdx[85] = 1.0*Reaction12_k4;
    dwdx[86] = 1.0*Reaction14_d5;
    dwdx[87] = 1.0*Reaction15_k5;
    dwdx[88] = 1.0*MEKPH*Reaction16_a6;
    dwdx[89] = 1.0*MAPK*Reaction19_a7;
    dwdx[90] = 1.0*MAPKp*Reaction25_a9;
    dwdx[91] = 1.0*Reaction17_d6;
    dwdx[92] = 1.0*Reaction18_k6;
    dwdx[93] = 1.0*RAFK*Reaction1_a1;
    dwdx[94] = 1.0*RAF*Reaction1_a1;
    dwdx[95] = 1.0*RAFp*Reaction4_a2;
    dwdx[96] = 1.0*Reaction2_d1;
    dwdx[97] = 1.0*Reaction3_k1;
    dwdx[98] = 1.0*MEKp*Reaction13_a5;
    dwdx[99] = 1.0*RAFPH*Reaction4_a2;
    dwdx[100] = 1.0*MEK*Reaction7_a3;
    dwdx[101] = 1.0*C2*kr1;
    dwdx[102] = 1.0*C2*kr1;
    dwdx[103] = 1.0*C6*kr1;
    dwdx[104] = 1.0*C6*kr1;
    dwdx[105] = 1.0*C9*kr1;
    dwdx[106] = 1.0*C9*kr1;
    dwdx[107] = 1.0*Reaction5_d2;
    dwdx[108] = 1.0*Reaction6_k2;
}