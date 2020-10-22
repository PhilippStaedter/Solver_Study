#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_Morris2002_CellCycle_CDK2Cyclin(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = CDK2cycA;
    x_solver[1] = CDK2cycA_star_;
    x_solver[2] = Cdk2;
    x_solver[3] = CyclinA;
}