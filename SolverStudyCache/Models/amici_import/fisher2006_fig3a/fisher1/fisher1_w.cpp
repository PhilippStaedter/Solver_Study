#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_fisher1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = -1.13e-13*NFAT_Nuc*k2 + 1.13e-13*NFAT_Pi_Nuc*k1;
    w[1] = 2.6900000000000001e-13*NFAT_Pi_Cyt*k3 - 1.13e-13*NFAT_Pi_Nuc*k4;
    w[2] = -2.6900000000000001e-13*Act_C_Cyt*NFAT_Cyt*k16 + 2.6900000000000001e-13*NFAT_Act_C_Cyt*k15;
    w[3] = -2.6900000000000001e-13*NFAT_Cyt*k2 + 2.6900000000000001e-13*NFAT_Pi_Cyt*k1;
    w[4] = -2.6900000000000001e-13*Act_C_Cyt*k20 + 2.6900000000000001e-13*pow(Ca_Cyt, 3)*Inact_C_Cyt*k19;
    w[5] = -1.13e-13*Act_C_Nuc*k20 + 1.13e-13*pow(Ca_Nuc, 3)*Inact_C_Nuc*k19;
    w[6] = 2.6900000000000001e-13*Inact_C_Cyt*k5 - 1.13e-13*Inact_C_Nuc*k6;
    w[7] = 2.6900000000000001e-13*Ca_Cyt*k21 - 1.13e-13*Ca_Nuc*k22;
    w[8] = 2.6900000000000001e-13*NFAT_Pi_Act_C_Cyt*k7 - 1.13e-13*NFAT_Pi_Act_C_Nuc*k8;
    w[9] = 1.13e-13*Act_C_Nuc*NFAT_Nuc*k16 - 1.13e-13*NFAT_Act_C_Nuc*k15;
    w[10] = -2.6900000000000001e-13*NFAT_Cyt*k17 + 1.13e-13*NFAT_Nuc*k18;
    w[11] = -2.6900000000000001e-13*Act_C_Cyt*k5 + 1.13e-13*Act_C_Nuc*k6;
    w[12] = 1.13e-13*NFAT_Act_C_Nuc*k14 - 1.13e-13*NFAT_Pi_Act_C_Nuc*k13;
    w[13] = -1.13e-13*Act_C_Nuc*NFAT_Pi_Nuc*k11 + 1.13e-13*NFAT_Pi_Act_C_Nuc*k12;
    w[14] = -2.6900000000000001e-13*NFAT_Act_C_Cyt*k9 + 1.13e-13*NFAT_Act_C_Nuc*k10;
    w[15] = 2.6900000000000001e-13*NFAT_Act_C_Cyt*k14 - 2.6900000000000001e-13*NFAT_Pi_Act_C_Cyt*k13;
    w[16] = -2.6900000000000001e-13*Act_C_Cyt*NFAT_Pi_Cyt*k11 + 2.6900000000000001e-13*NFAT_Pi_Act_C_Cyt*k12;
}