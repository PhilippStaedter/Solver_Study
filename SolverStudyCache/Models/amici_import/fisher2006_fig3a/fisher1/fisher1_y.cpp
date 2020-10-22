#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_fisher1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Act_C_Cyt;
    y[1] = Act_C_Nuc;
    y[2] = Ca_Cyt;
    y[3] = Ca_Nuc;
    y[4] = Inact_C_Cyt;
    y[5] = Inact_C_Nuc;
    y[6] = NFAT_Act_C_Cyt;
    y[7] = NFAT_Act_C_Nuc;
    y[8] = NFAT_Cyt;
    y[9] = NFAT_Nuc;
    y[10] = NFAT_Pi_Act_C_Cyt;
    y[11] = NFAT_Pi_Act_C_Nuc;
    y[12] = NFAT_Pi_Cyt;
    y[13] = NFAT_Pi_Nuc;
}