#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_kofahl(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Bar1;
    y[1] = Bar1a;
    y[2] = Bar1aex;
    y[3] = Cdc28;
    y[4] = Far1;
    y[5] = Far1PP;
    y[6] = Far1U;
    y[7] = Fus3;
    y[8] = Fus3PP;
    y[9] = GaGDP;
    y[10] = GaGTP;
    y[11] = Gabc;
    y[12] = Gbc;
    y[13] = Sst2;
    y[14] = Ste11;
    y[15] = Ste12;
    y[16] = Ste12a;
    y[17] = Ste2;
    y[18] = Ste20;
    y[19] = Ste2a;
    y[20] = Ste5;
    y[21] = Ste7;
    y[22] = alpha;
    y[23] = complexA;
    y[24] = complexB;
    y[25] = complexC;
    y[26] = complexD;
    y[27] = complexE;
    y[28] = complexF;
    y[29] = complexG;
    y[30] = complexH;
    y[31] = complexI;
    y[32] = complexK;
    y[33] = complexL;
    y[34] = complexM;
    y[35] = complexN;
    y[36] = amici_p;
}