#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_hockin2(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = APC;
    x_solver[1] = HC;
    x_solver[2] = HC3;
    x_solver[3] = HC36;
    x_solver[4] = HC5;
    x_solver[5] = HC53;
    x_solver[6] = HC536;
    x_solver[7] = HC56;
    x_solver[8] = LC;
    x_solver[9] = LC_APC;
    x_solver[10] = Va;
    x_solver[11] = Va3;
    x_solver[12] = Va36;
    x_solver[13] = Va36_APC;
    x_solver[14] = Va3_APC;
    x_solver[15] = Va5;
    x_solver[16] = Va53;
    x_solver[17] = Va536;
    x_solver[18] = Va536_APC;
    x_solver[19] = Va53_APC;
    x_solver[20] = Va56;
    x_solver[21] = Va56_APC;
    x_solver[22] = Va5_APC;
    x_solver[23] = VaA3;
    x_solver[24] = VaA36;
    x_solver[25] = VaA53;
    x_solver[26] = VaA536;
    x_solver[27] = VaLCA1;
    x_solver[28] = VaLCA1_APC;
    x_solver[29] = Va_APC;
}