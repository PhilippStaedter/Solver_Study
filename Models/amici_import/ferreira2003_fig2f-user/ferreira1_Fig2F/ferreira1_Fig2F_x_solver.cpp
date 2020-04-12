#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_ferreira1_Fig2F(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Amadori;
    x_solver[1] = CML;
    x_solver[2] = Glucose;
    x_solver[3] = Glyoxal;
    x_solver[4] = Lysine;
    x_solver[5] = Schiff;
}