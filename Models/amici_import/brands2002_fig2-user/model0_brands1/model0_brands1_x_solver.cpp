#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_brands1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = AMP;
    x_solver[1] = Acetic_acid;
    x_solver[2] = Amadori;
    x_solver[3] = C5;
    x_solver[4] = Cn;
    x_solver[5] = Formic_acid;
    x_solver[6] = Fru;
    x_solver[7] = Glu;
    x_solver[8] = Melanoidin;
    x_solver[9] = Triose;
    x_solver[10] = lys_R;
}