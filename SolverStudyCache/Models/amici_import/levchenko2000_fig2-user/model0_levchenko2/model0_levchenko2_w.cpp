#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_levchenko2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*RAF*RAFK*Reaction1_a1;
    w[1] = 1.0*MEKPH*MEKp*Reaction10_a4;
    w[2] = 1.0*MEKpMEKPH*Reaction11_d4;
    w[3] = 1.0*MEKpMEKPH*Reaction12_k4;
    w[4] = 1.0*MEKp*RAFp*Reaction13_a5;
    w[5] = 1.0*MEKpRAFp*Reaction14_d5;
    w[6] = 1.0*MEKpRAFp*Reaction15_k5;
    w[7] = 1.0*MEKPH*MEKpp*Reaction16_a6;
    w[8] = 1.0*MEKppMEKPH*Reaction17_d6;
    w[9] = 1.0*MEKppMEKPH*Reaction18_k6;
    w[10] = 1.0*MAPK*MEKpp*Reaction19_a7;
    w[11] = 1.0*RAFRAFK*Reaction2_d1;
    w[12] = 1.0*MAPKMEKpp*Reaction20_d7;
    w[13] = 1.0*MAPKMEKpp*Reaction21_k7;
    w[14] = 1.0*MAPKPH*MAPKp*Reaction22_a8;
    w[15] = 1.0*MAPKpMAPKPH*Reaction23_d8;
    w[16] = 1.0*MAPKpMAPKPH*Reaction24_k8;
    w[17] = 1.0*MAPKp*MEKpp*Reaction25_a9;
    w[18] = 1.0*MAPKpMEKpp*Reaction26_d9;
    w[19] = 1.0*MAPKpMEKpp*Reaction27_k9;
    w[20] = 1.0*MAPKPH*MAPKpp*Reaction28_a10;
    w[21] = 1.0*MAPKppMAPKPH*Reaction29_d10;
    w[22] = 1.0*RAFRAFK*Reaction3_k1;
    w[23] = 1.0*MAPKppMAPKPH*Reaction30_k10;
    w[24] = 1.0*RAFPH*RAFp*Reaction4_a2;
    w[25] = 1.0*RAFpRAFPH*Reaction5_d2;
    w[26] = 1.0*RAFpRAFPH*Reaction6_k2;
    w[27] = 1.0*MEK*RAFp*Reaction7_a3;
    w[28] = 1.0*MEKRAFp*Reaction8_d3;
    w[29] = 1.0*MEKRAFp*Reaction9_k3;
    w[30] = C1*(1.0*MAPK*on2 + 1.0*MEK*on1);
    w[31] = C2*of1 + C3*of3 + C4*of2 + C5*of4;
    w[32] = 1.0*C2*RAFp*kr1 + C2*(1.0*MAPK*on2 + of1);
    w[33] = 1.0*C1*MEK*on1 + C6*of2 + C9*of4;
    w[34] = 1.0*C3*MAPK*on2 + C3*of3;
    w[35] = 1.0*C2*RAFp*kr1 + C7*of2 + C8*of4;
    w[36] = C4*(1.0*MEK*on1 + of2);
    w[37] = 1.0*C1*MAPK*on2 + C6*of1 + C7*of3;
    w[38] = C5*(1.0*MEK*on1 + of4);
    w[39] = C8*of3 + C9*of1;
    w[40] = 1.0*C6*RAFp*kr1 + C6*(of1 + of2);
    w[41] = 1.0*C2*MAPK*on2 + 1.0*C4*MEK*on1;
    w[42] = C7*kr2 + C7*of2 + C7*of3;
    w[43] = 1.0*C3*MAPK*on2 + 1.0*C6*RAFp*kr1;
    w[44] = C8*(of3 + of4);
    w[45] = C7*kr2 + 1.0*C9*RAFp*kr1;
    w[46] = 1.0*C9*RAFp*kr1 + C9*(of1 + of4);
    w[47] = 1.0*C5*MEK*on1;
    w[48] = of4*(C5 + C8 + C9);
    w[49] = 1.0*MAPK*on2*(C1 + C2 + C3);
    w[50] = of2*(C4 + C6 + C7);
    w[51] = of3*(C3 + C7 + C8);
    w[52] = 1.0*MEK*on1*(C1 + C4 + C5);
    w[53] = of1*(C2 + C6 + C9);
}