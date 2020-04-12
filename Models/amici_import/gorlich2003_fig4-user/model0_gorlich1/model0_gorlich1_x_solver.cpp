#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_gorlich1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = GDP;
    x_solver[1] = GTP;
    x_solver[2] = RCC1;
    x_solver[3] = RCC1_Ran;
    x_solver[4] = RCC1_RanGDP;
    x_solver[5] = RCC1_RanGTP;
    x_solver[6] = RanBP1;
    x_solver[7] = RanGAP;
    x_solver[8] = RanGDP_cy;
    x_solver[9] = RanGDP_nuc;
    x_solver[10] = RanGTP_RanBP1;
    x_solver[11] = RanGTP_cy;
    x_solver[12] = RanGTP_nuc;
}