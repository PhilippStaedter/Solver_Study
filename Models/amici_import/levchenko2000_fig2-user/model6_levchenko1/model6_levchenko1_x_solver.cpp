#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model6_levchenko1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = MAPK;
    x_solver[1] = MAPKMEKpp;
    x_solver[2] = MAPKPH;
    x_solver[3] = MAPKp;
    x_solver[4] = MAPKpMAPKPH;
    x_solver[5] = MAPKpMEKpp;
    x_solver[6] = MAPKpp;
    x_solver[7] = MAPKppMAPKPH;
    x_solver[8] = MEK;
    x_solver[9] = MEKPH;
    x_solver[10] = MEKRAFp;
    x_solver[11] = MEKp;
    x_solver[12] = MEKpMEKPH;
    x_solver[13] = MEKpRAFp;
    x_solver[14] = MEKpp;
    x_solver[15] = MEKppMEKPH;
    x_solver[16] = RAF;
    x_solver[17] = RAFK;
    x_solver[18] = RAFPH;
    x_solver[19] = RAFRAFK;
    x_solver[20] = RAFp;
    x_solver[21] = RAFpRAFPH;
}