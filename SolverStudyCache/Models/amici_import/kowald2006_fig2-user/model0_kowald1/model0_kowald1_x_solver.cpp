#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_kowald1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = species_0000001;
    x_solver[1] = species_0000002;
    x_solver[2] = species_0000006;
    x_solver[3] = species_0000007;
    x_solver[4] = species_0000008;
    x_solver[5] = species_0000009;
    x_solver[6] = species_0000011;
    x_solver[7] = species_0000016;
    x_solver[8] = species_0000017;
}