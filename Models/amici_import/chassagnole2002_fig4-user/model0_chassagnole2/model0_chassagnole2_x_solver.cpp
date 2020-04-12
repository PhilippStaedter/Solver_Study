#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_chassagnole2(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = cdhap;
    x_solver[1] = ce4p;
    x_solver[2] = cf6p;
    x_solver[3] = cfdp;
    x_solver[4] = cg1p;
    x_solver[5] = cg6p;
    x_solver[6] = cgap;
    x_solver[7] = cglcex;
    x_solver[8] = cpep;
    x_solver[9] = cpg;
    x_solver[10] = cpg2;
    x_solver[11] = cpg3;
    x_solver[12] = cpgp;
    x_solver[13] = cpyr;
    x_solver[14] = crib5p;
    x_solver[15] = cribu5p;
    x_solver[16] = csed7p;
    x_solver[17] = cxyl5p;
}