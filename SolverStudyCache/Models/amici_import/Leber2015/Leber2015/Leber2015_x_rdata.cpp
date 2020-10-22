#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_Leber2015(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Cdiff;
    x_rdata[1] = Commensal_Beneficial;
    x_rdata[2] = Commensal_Dead;
    x_rdata[3] = tDC_LP;
    x_rdata[4] = tDC_MLN;
    x_rdata[5] = Commensal_Harmful;
    x_rdata[6] = N_Lum;
    x_rdata[7] = E;
    x_rdata[8] = E_d;
    x_rdata[9] = iDC_E;
    x_rdata[10] = E_i;
    x_rdata[11] = M_LP;
    x_rdata[12] = eDC_LP;
    x_rdata[13] = M0;
    x_rdata[14] = N_LP;
    x_rdata[15] = Th17_LP;
    x_rdata[16] = Th1_LP;
    x_rdata[17] = iTreg_LP;
    x_rdata[18] = eDC_MLN;
    x_rdata[19] = iTreg_MLN;
    x_rdata[20] = nT;
    x_rdata[21] = Th17_MLN;
    x_rdata[22] = Th1_MLN;
}