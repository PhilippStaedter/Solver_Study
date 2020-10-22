#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_kofahl(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Bar1;
    x_solver[1] = Bar1a;
    x_solver[2] = Bar1aex;
    x_solver[3] = Cdc28;
    x_solver[4] = Far1;
    x_solver[5] = Far1PP;
    x_solver[6] = Far1U;
    x_solver[7] = Fus3;
    x_solver[8] = Fus3PP;
    x_solver[9] = GaGDP;
    x_solver[10] = GaGTP;
    x_solver[11] = Gabc;
    x_solver[12] = Gbc;
    x_solver[13] = Sst2;
    x_solver[14] = Ste11;
    x_solver[15] = Ste12;
    x_solver[16] = Ste12a;
    x_solver[17] = Ste2;
    x_solver[18] = Ste20;
    x_solver[19] = Ste2a;
    x_solver[20] = Ste5;
    x_solver[21] = Ste7;
    x_solver[22] = alpha;
    x_solver[23] = complexA;
    x_solver[24] = complexB;
    x_solver[25] = complexC;
    x_solver[26] = complexD;
    x_solver[27] = complexE;
    x_solver[28] = complexF;
    x_solver[29] = complexG;
    x_solver[30] = complexH;
    x_solver[31] = complexI;
    x_solver[32] = complexK;
    x_solver[33] = complexL;
    x_solver[34] = complexM;
    x_solver[35] = complexN;
    x_solver[36] = amici_p;
}