#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_kofahl(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = Ste12a*k36;
    dwdx[1] = k37;
    dwdx[2] = k38;
    dwdx[3] = alpha*k1;
    dwdx[4] = Far1*k41;
    dwdx[5] = Far1PP*k45;
    dwdx[6] = pow(Fus3PP, 2)*k39/(pow(Fus3PP, 2) + 10000);
    dwdx[7] = Cdc28*k41;
    dwdx[8] = k40;
    dwdx[9] = Gbc*k42;
    dwdx[10] = Cdc28*k45;
    dwdx[11] = Ste7*k14;
    dwdx[12] = complexL*k29;
    dwdx[13] = k33;
    dwdx[14] = Ste12*k34;
    dwdx[15] = -2*Far1*pow(Fus3PP, 3)*k39/pow(pow(Fus3PP, 2) + 10000, 2) + 2*Far1*Fus3PP*k39/(pow(Fus3PP, 2) + 10000);
    dwdx[16] = -2*pow(Fus3PP, 3)*k46/pow(pow(Fus3PP, 2) + 16, 2) + 2*Fus3PP*k46/(pow(Fus3PP, 2) + 16);
    dwdx[17] = Gbc*k9;
    dwdx[18] = k7;
    dwdx[19] = Sst2*k8;
    dwdx[20] = Ste2a*k6;
    dwdx[21] = complexC*k10;
    dwdx[22] = Far1PP*k42;
    dwdx[23] = GaGDP*k9;
    dwdx[24] = k47;
    dwdx[25] = GaGTP*k8;
    dwdx[26] = Ste5*k12;
    dwdx[27] = Fus3PP*k34;
    dwdx[28] = k35;
    dwdx[29] = Bar1*k36;
    dwdx[30] = alpha*k2;
    dwdx[31] = k5;
    dwdx[32] = complexD*k18;
    dwdx[33] = k3;
    dwdx[34] = k4;
    dwdx[35] = Gabc*k6;
    dwdx[36] = Ste11*k12;
    dwdx[37] = Fus3*k14;
    dwdx[38] = Bar1aex*k1;
    dwdx[39] = Ste2*k2;
    dwdx[40] = k13;
    dwdx[41] = complexB*k16;
    dwdx[42] = k15;
    dwdx[43] = complexA*k16;
    dwdx[44] = Gbc*k10;
    dwdx[45] = k17;
    dwdx[46] = k11;
    dwdx[47] = Ste20*k18;
    dwdx[48] = k19;
    dwdx[49] = k20;
    dwdx[50] = k21;
    dwdx[51] = k22;
    dwdx[52] = k23;
    dwdx[53] = k24;
    dwdx[54] = k25;
    dwdx[55] = k26;
    dwdx[56] = k27;
    dwdx[57] = k28;
    dwdx[58] = k30;
    dwdx[59] = k31;
    dwdx[60] = Fus3*k29;
    dwdx[61] = k32;
    dwdx[62] = k43;
    dwdx[63] = k44;
}