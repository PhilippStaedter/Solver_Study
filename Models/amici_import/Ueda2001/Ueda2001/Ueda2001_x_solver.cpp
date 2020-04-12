#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_Ueda2001(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = EmptySet;
    x_solver[1] = CCc;
    x_solver[2] = CCn;
    x_solver[3] = Clkc;
    x_solver[4] = Clkm;
    x_solver[5] = Perc;
    x_solver[6] = Perm;
    x_solver[7] = PTc;
    x_solver[8] = PTn;
    x_solver[9] = Timc;
    x_solver[10] = Timm;
    x_solver[11] = species_0000012;
    x_solver[12] = species_0000013;
}