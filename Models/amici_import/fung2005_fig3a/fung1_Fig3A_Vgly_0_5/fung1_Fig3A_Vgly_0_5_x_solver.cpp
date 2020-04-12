#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_fung1_Fig3A_Vgly_0_5(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = AcCoA;
    x_solver[1] = AcP;
    x_solver[2] = Acs;
    x_solver[3] = HOAc;
    x_solver[4] = HOAc_E;
    x_solver[5] = LacI;
    x_solver[6] = OAc;
    x_solver[7] = Pta;
}