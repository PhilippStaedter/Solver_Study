#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_hatakeyama1_Fig5G(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Akt;
    x_solver[1] = AktPIP;
    x_solver[2] = AktPIP3;
    x_solver[3] = AktPIPP;
    x_solver[4] = E;
    x_solver[5] = ERK;
    x_solver[6] = ERKP;
    x_solver[7] = ERKPP;
    x_solver[8] = GS;
    x_solver[9] = HRG;
    x_solver[10] = MEK;
    x_solver[11] = MEKP;
    x_solver[12] = MEKPP;
    x_solver[13] = MKP3;
    x_solver[14] = PI3K;
    x_solver[15] = PI3Kstar;
    x_solver[16] = PIP3;
    x_solver[17] = PP2A;
    x_solver[18] = P_I;
    x_solver[19] = R;
    x_solver[20] = RHRG;
    x_solver[21] = RHRG2;
    x_solver[22] = RP;
    x_solver[23] = RPI3K;
    x_solver[24] = RPI3Kstar;
    x_solver[25] = RShGS;
    x_solver[26] = RShP;
    x_solver[27] = RShc;
    x_solver[28] = Raf;
    x_solver[29] = Rafstar;
    x_solver[30] = RasGDP;
    x_solver[31] = RasGTP;
    x_solver[32] = ShGS;
    x_solver[33] = ShP;
    x_solver[34] = Shc;
    x_solver[35] = internalization;
}