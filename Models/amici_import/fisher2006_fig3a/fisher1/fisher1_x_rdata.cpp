#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_fisher1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Act_C_Cyt;
    x_rdata[1] = Act_C_Nuc;
    x_rdata[2] = Ca_Cyt;
    x_rdata[3] = Ca_Nuc;
    x_rdata[4] = Inact_C_Cyt;
    x_rdata[5] = Inact_C_Nuc;
    x_rdata[6] = NFAT_Act_C_Cyt;
    x_rdata[7] = NFAT_Act_C_Nuc;
    x_rdata[8] = NFAT_Cyt;
    x_rdata[9] = NFAT_Nuc;
    x_rdata[10] = NFAT_Pi_Act_C_Cyt;
    x_rdata[11] = NFAT_Pi_Act_C_Nuc;
    x_rdata[12] = NFAT_Pi_Cyt;
    x_rdata[13] = NFAT_Pi_Nuc;
}