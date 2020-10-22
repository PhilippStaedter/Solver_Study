#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_aguda1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = CYCDcdk4;
    y[1] = CYCDcdk4p27;
    y[2] = CYCEcdk2p27;
    y[3] = E2F;
    y[4] = aCYCEcdk2;
    y[5] = aCYCEcdk20;
    y[6] = aCYCEcdk21;
    y[7] = iCYCEcdk2;
    y[8] = p27;
    y[9] = pRB;
    y[10] = pRBE2F;
    y[11] = xvar;
}