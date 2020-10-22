#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model2_mcauley1(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[1] = -1.0*dwdx1;
    JDiag[2] = -1.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6 + 1.0*dwdx7;
    JDiag[4] = -1.0*dwdx9;
    JDiag[8] = -1.0*dwdx14 - 1.0*dwdx15;
    JDiag[9] = -1.0*dwdx16;
    JDiag[11] = 1.0*dwdx19 - 1.0*dwdx20 - 1.0*dwdx21;
    JDiag[13] = -1.0*dwdx22 - 1.0*dwdx23;
    JDiag[15] = -1.0*dwdx25 - 1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28;
    JDiag[17] = -1.0*dwdx31;
    JDiag[20] = -1.0*dwdx33;
    JDiag[23] = -1.0*dwdx34 - 1.0*dwdx35 - 1.0*dwdx36;
    JDiag[28] = -1.0*dwdx41 + 1.0*dwdx42;
    JDiag[29] = -1.0*dwdx43 - 1.0*dwdx44;
    JDiag[31] = -1.0*dwdx47 + 1.0*dwdx48 - 1.0*dwdx49 - 1.0*dwdx50 - 1.0*dwdx52;
}