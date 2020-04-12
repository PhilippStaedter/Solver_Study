#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model0_sarma3(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0 + 1.0*dwdx1;
    JDiag[3] = -1.0*dwdx9;
    JDiag[4] = -1.0*dwdx16 + 1.0*dwdx19;
    JDiag[5] = 1.0*dwdx23 - 1.0*dwdx24 + 1.0*dwdx25 - 1.0*dwdx26;
    JDiag[6] = -1.0*dwdx30;
    JDiag[7] = 1.0*dwdx34 - 1.0*dwdx37;
    JDiag[8] = -1.0*dwdx40 + 1.0*dwdx43 - 1.0*dwdx44 + 1.0*dwdx45;
    JDiag[9] = -1.0*dwdx49;
}