#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_levchenko2(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = C1;
    x_solver[1] = C2;
    x_solver[2] = C3;
    x_solver[3] = C4;
    x_solver[4] = C5;
    x_solver[5] = C6;
    x_solver[6] = C7;
    x_solver[7] = C8;
    x_solver[8] = C9;
    x_solver[9] = MAPK;
    x_solver[10] = MAPKMEKpp;
    x_solver[11] = MAPKPH;
    x_solver[12] = MAPKp;
    x_solver[13] = MAPKpMAPKPH;
    x_solver[14] = MAPKpMEKpp;
    x_solver[15] = MAPKpp;
    x_solver[16] = MAPKppMAPKPH;
    x_solver[17] = MEK;
    x_solver[18] = MEKPH;
    x_solver[19] = MEKRAFp;
    x_solver[20] = MEKp;
    x_solver[21] = MEKpMEKPH;
    x_solver[22] = MEKpRAFp;
    x_solver[23] = MEKpp;
    x_solver[24] = MEKppMEKPH;
    x_solver[25] = RAF;
    x_solver[26] = RAFK;
    x_solver[27] = RAFPH;
    x_solver[28] = RAFRAFK;
    x_solver[29] = RAFp;
    x_solver[30] = RAFpRAFPH;
}