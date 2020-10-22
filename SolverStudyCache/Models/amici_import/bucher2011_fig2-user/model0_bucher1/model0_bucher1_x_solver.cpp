#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_bucher1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = ASL_b;
    x_solver[1] = ASL_c;
    x_solver[2] = ASL_m;
    x_solver[3] = ASLoOH_b;
    x_solver[4] = ASLoOH_c;
    x_solver[5] = ASLoOH_m;
    x_solver[6] = ASLpOH_b;
    x_solver[7] = ASLpOH_c;
    x_solver[8] = ASLpOH_m;
    x_solver[9] = AS_b;
    x_solver[10] = AS_c;
    x_solver[11] = AS_m;
    x_solver[12] = ASoOH_b;
    x_solver[13] = ASoOH_c;
    x_solver[14] = ASoOH_m;
    x_solver[15] = ASpOH_b;
    x_solver[16] = ASpOH_c;
    x_solver[17] = ASpOH_m;
}