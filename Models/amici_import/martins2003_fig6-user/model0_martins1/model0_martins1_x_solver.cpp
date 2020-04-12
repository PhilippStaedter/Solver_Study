#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_martins1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = AA;
    x_solver[1] = Cn;
    x_solver[2] = DFG;
    x_solver[3] = E1;
    x_solver[4] = E2;
    x_solver[5] = FA;
    x_solver[6] = Fru;
    x_solver[7] = Glu;
    x_solver[8] = Gly;
    x_solver[9] = MG;
    x_solver[10] = Man;
    x_solver[11] = Mel;
    x_solver[12] = _1DG;
    x_solver[13] = _3DG;
}