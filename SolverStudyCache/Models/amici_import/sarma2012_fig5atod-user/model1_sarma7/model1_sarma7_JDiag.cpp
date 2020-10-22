#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model1_sarma7(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0;
    JDiag[1] = -1.0*dwdx3;
    JDiag[3] = -1.0*dwdx6;
    JDiag[4] = -1.0*dwdx7;
    JDiag[5] = -1.0*dwdx8;
    JDiag[6] = 1.0*dwdx10 - 1.0*dwdx11 + 1.0*dwdx12 - 1.0*dwdx13;
    JDiag[7] = -1.0*dwdx16;
    JDiag[8] = -1.0*dwdx18;
    JDiag[9] = 1.0*dwdx20 - 1.0*dwdx21 + 1.0*dwdx22 - 1.0*dwdx23;
    JDiag[10] = -1.0*dwdx28;
}