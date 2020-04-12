#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model1_mellor1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = species_1;
    x_solver[1] = species_10;
    x_solver[2] = species_11;
    x_solver[3] = species_12;
    x_solver[4] = species_13;
    x_solver[5] = species_14;
    x_solver[6] = species_15;
    x_solver[7] = species_7;
    x_solver[8] = species_8;
    x_solver[9] = species_9;
}