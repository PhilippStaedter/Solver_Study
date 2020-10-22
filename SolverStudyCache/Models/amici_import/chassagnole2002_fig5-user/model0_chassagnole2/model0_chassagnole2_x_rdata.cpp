#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_chassagnole2(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = cdhap;
    x_rdata[1] = ce4p;
    x_rdata[2] = cf6p;
    x_rdata[3] = cfdp;
    x_rdata[4] = cg1p;
    x_rdata[5] = cg6p;
    x_rdata[6] = cgap;
    x_rdata[7] = cglcex;
    x_rdata[8] = cpep;
    x_rdata[9] = cpg;
    x_rdata[10] = cpg2;
    x_rdata[11] = cpg3;
    x_rdata[12] = cpgp;
    x_rdata[13] = cpyr;
    x_rdata[14] = crib5p;
    x_rdata[15] = cribu5p;
    x_rdata[16] = csed7p;
    x_rdata[17] = cxyl5p;
}