#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_kholodenko1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = EGF;
    x_solver[1] = GS;
    x_solver[2] = Grb;
    x_solver[3] = PLCg;
    x_solver[4] = PLCgP;
    x_solver[5] = PLCgl;
    x_solver[6] = R;
    x_solver[7] = R2;
    x_solver[8] = RG;
    x_solver[9] = RGS;
    x_solver[10] = RP;
    x_solver[11] = RPLCg;
    x_solver[12] = RPLCgP;
    x_solver[13] = RSh;
    x_solver[14] = RShG;
    x_solver[15] = RShGS;
    x_solver[16] = RShP;
    x_solver[17] = Ra;
    x_solver[18] = SOS;
    x_solver[19] = ShG;
    x_solver[20] = ShGS;
    x_solver[21] = ShP;
    x_solver[22] = Shc;
}