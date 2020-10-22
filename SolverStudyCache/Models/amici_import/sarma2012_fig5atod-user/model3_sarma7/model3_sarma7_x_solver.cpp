#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model3_sarma7(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = species_0;
    x_solver[1] = species_1;
    x_solver[2] = species_10;
    x_solver[3] = species_11;
    x_solver[4] = species_12;
    x_solver[5] = species_2;
    x_solver[6] = species_3;
    x_solver[7] = species_4;
    x_solver[8] = species_5;
    x_solver[9] = species_6;
    x_solver[10] = species_7;
    x_solver[11] = species_8;
    x_solver[12] = species_9;
}