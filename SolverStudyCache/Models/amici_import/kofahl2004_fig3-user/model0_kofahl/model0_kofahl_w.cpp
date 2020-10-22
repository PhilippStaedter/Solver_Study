#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_kofahl(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = Bar1aex*alpha*k1;
    w[1] = Gbc*complexC*k10;
    w[2] = complexD*k11;
    w[3] = Ste11*Ste5*k12;
    w[4] = complexA*k13;
    w[5] = Fus3*Ste7*k14;
    w[6] = complexB*k15;
    w[7] = complexA*complexB*k16;
    w[8] = complexC*k17;
    w[9] = Ste20*complexD*k18;
    w[10] = complexE*k19;
    w[11] = Ste2*alpha*k2;
    w[12] = complexE*k20;
    w[13] = complexE*k21;
    w[14] = complexF*k22;
    w[15] = complexF*k23;
    w[16] = complexG*k24;
    w[17] = complexG*k25;
    w[18] = complexH*k26;
    w[19] = complexH*k27;
    w[20] = complexI*k28;
    w[21] = Fus3*complexL*k29;
    w[22] = Ste2a*k3;
    w[23] = complexK*k30;
    w[24] = complexK*k31;
    w[25] = complexL*k32;
    w[26] = Fus3PP*k33;
    w[27] = Fus3PP*Ste12*k34;
    w[28] = Ste12a*k35;
    w[29] = Bar1*Ste12a*k36;
    w[30] = Bar1a*k37;
    w[31] = Bar1a*k38;
    w[32] = Far1*pow(Fus3PP, 2)*k39/(pow(Fus3PP, 2) + 10000);
    w[33] = Ste2a*k4;
    w[34] = Far1PP*k40;
    w[35] = Cdc28*Far1*k41;
    w[36] = Far1PP*Gbc*k42;
    w[37] = complexM*k43;
    w[38] = complexN*k44;
    w[39] = Cdc28*Far1PP*k45;
    w[40] = pow(Fus3PP, 2)*k46/(pow(Fus3PP, 2) + 16);
    w[41] = Sst2*k47;
    w[42] = Ste2*k5;
    w[43] = Gabc*Ste2a*k6;
    w[44] = GaGTP*k7;
    w[45] = GaGTP*Sst2*k8;
    w[46] = GaGDP*Gbc*k9;
}