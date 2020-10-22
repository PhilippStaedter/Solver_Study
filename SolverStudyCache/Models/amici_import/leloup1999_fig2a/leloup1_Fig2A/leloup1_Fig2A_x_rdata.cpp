#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_leloup1_Fig2A(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = CC;
    x_rdata[1] = Cn;
    x_rdata[2] = Mp;
    x_rdata[3] = Mt;
    x_rdata[4] = P0;
    x_rdata[5] = P1;
    x_rdata[6] = P2;
    x_rdata[7] = T0;
    x_rdata[8] = T1;
    x_rdata[9] = T2;
}