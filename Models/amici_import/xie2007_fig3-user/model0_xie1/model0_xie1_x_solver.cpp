#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_xie1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = CC;
    x_solver[1] = CCPT;
    x_solver[2] = CLK;
    x_solver[3] = CYC;
    x_solver[4] = PDP;
    x_solver[5] = PER;
    x_solver[6] = PT;
    x_solver[7] = TIM;
    x_solver[8] = VRI;
    x_solver[9] = clkm;
    x_solver[10] = clkp;
    x_solver[11] = pdpm;
    x_solver[12] = pdpp;
    x_solver[13] = perm;
    x_solver[14] = perp;
    x_solver[15] = prcpdp;
    x_solver[16] = prcper;
    x_solver[17] = prct;
    x_solver[18] = prcv;
    x_solver[19] = prpc;
    x_solver[20] = prvc;
    x_solver[21] = timm;
    x_solver[22] = timp;
    x_solver[23] = vrim;
    x_solver[24] = vrip;
}