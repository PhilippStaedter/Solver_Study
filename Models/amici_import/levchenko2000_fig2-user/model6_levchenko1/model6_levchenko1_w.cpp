#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model6_levchenko1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
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
}