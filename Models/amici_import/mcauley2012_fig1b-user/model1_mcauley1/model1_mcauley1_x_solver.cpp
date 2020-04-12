#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model1_mcauley1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = species_1;
    x_solver[1] = species_10;
    x_solver[2] = species_11;
    x_solver[3] = species_12;
    x_solver[4] = species_13;
    x_solver[5] = species_14;
    x_solver[6] = species_15;
    x_solver[7] = species_16;
    x_solver[8] = species_17;
    x_solver[9] = species_18;
    x_solver[10] = species_19;
    x_solver[11] = species_2;
    x_solver[12] = species_20;
    x_solver[13] = species_21;
    x_solver[14] = species_22;
    x_solver[15] = species_23;
    x_solver[16] = species_24;
    x_solver[17] = species_25;
    x_solver[18] = species_26;
    x_solver[19] = species_27;
    x_solver[20] = species_28;
    x_solver[21] = species_29;
    x_solver[22] = species_3;
    x_solver[23] = species_30;
    x_solver[24] = species_31;
    x_solver[25] = species_32;
    x_solver[26] = species_33;
    x_solver[27] = species_34;
    x_solver[28] = species_4;
    x_solver[29] = species_5;
    x_solver[30] = species_6;
    x_solver[31] = species_7;
    x_solver[32] = species_8;
    x_solver[33] = species_9;
}