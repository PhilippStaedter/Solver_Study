#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model2_(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = EtOH;
    x_solver[1] = Glycerol;
    x_solver[2] = PiVac;
    x_solver[3] = atp;
    x_solver[4] = fbp;
    x_solver[5] = g6p;
    x_solver[6] = phos;
}