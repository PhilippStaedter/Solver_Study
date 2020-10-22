#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model6_levchenko1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*MEKpp*Reaction19_a7;
    dwdx[1] = 1.0*Reaction20_d7;
    dwdx[2] = 1.0*Reaction21_k7;
    dwdx[3] = 1.0*MAPKp*Reaction22_a8;
    dwdx[4] = 1.0*MAPKpp*Reaction28_a10;
    dwdx[5] = 1.0*MAPKPH*Reaction22_a8;
    dwdx[6] = 1.0*MEKpp*Reaction25_a9;
    dwdx[7] = 1.0*Reaction23_d8;
    dwdx[8] = 1.0*Reaction24_k8;
    dwdx[9] = 1.0*Reaction26_d9;
    dwdx[10] = 1.0*Reaction27_k9;
    dwdx[11] = 1.0*MAPKPH*Reaction28_a10;
    dwdx[12] = 1.0*Reaction29_d10;
    dwdx[13] = 1.0*Reaction30_k10;
    dwdx[14] = 1.0*RAFp*Reaction7_a3;
    dwdx[15] = 1.0*MEKp*Reaction10_a4;
    dwdx[16] = 1.0*MEKpp*Reaction16_a6;
    dwdx[17] = 1.0*Reaction8_d3;
    dwdx[18] = 1.0*Reaction9_k3;
    dwdx[19] = 1.0*MEKPH*Reaction10_a4;
    dwdx[20] = 1.0*RAFp*Reaction13_a5;
    dwdx[21] = 1.0*Reaction11_d4;
    dwdx[22] = 1.0*Reaction12_k4;
    dwdx[23] = 1.0*Reaction14_d5;
    dwdx[24] = 1.0*Reaction15_k5;
    dwdx[25] = 1.0*MEKPH*Reaction16_a6;
    dwdx[26] = 1.0*MAPK*Reaction19_a7;
    dwdx[27] = 1.0*MAPKp*Reaction25_a9;
    dwdx[28] = 1.0*Reaction17_d6;
    dwdx[29] = 1.0*Reaction18_k6;
    dwdx[30] = 1.0*RAFK*Reaction1_a1;
    dwdx[31] = 1.0*RAF*Reaction1_a1;
    dwdx[32] = 1.0*RAFp*Reaction4_a2;
    dwdx[33] = 1.0*Reaction2_d1;
    dwdx[34] = 1.0*Reaction3_k1;
    dwdx[35] = 1.0*MEKp*Reaction13_a5;
    dwdx[36] = 1.0*RAFPH*Reaction4_a2;
    dwdx[37] = 1.0*MEK*Reaction7_a3;
    dwdx[38] = 1.0*Reaction5_d2;
    dwdx[39] = 1.0*Reaction6_k2;
}