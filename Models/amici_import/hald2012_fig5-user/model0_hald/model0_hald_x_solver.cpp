#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_hald(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = ACA;
    x_solver[1] = ACA0;
    x_solver[2] = ADP;
    x_solver[3] = AMP;
    x_solver[4] = ATP;
    x_solver[5] = DHAP;
    x_solver[6] = DHAPCN;
    x_solver[7] = DPG;
    x_solver[8] = EtOH;
    x_solver[9] = EtOH0;
    x_solver[10] = F6P;
    x_solver[11] = FBP;
    x_solver[12] = G6P;
    x_solver[13] = GAP;
    x_solver[14] = Glc;
    x_solver[15] = Glc0;
    x_solver[16] = Glyc;
    x_solver[17] = Glyc0;
    x_solver[18] = HCN;
    x_solver[19] = HCN0;
    x_solver[20] = NAD;
    x_solver[21] = NADH;
    x_solver[22] = OAc;
    x_solver[23] = OAc0;
    x_solver[24] = PEP;
    x_solver[25] = Pyr;
    x_solver[26] = PyrCN;
    x_solver[27] = X;
    x_solver[28] = drain;
    x_solver[29] = glycogen;
    x_solver[30] = lacto;
}