#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model0_zi1(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -952.38095238095241*dwdx0 - 952.38095238095241*dwdx1;
    JDiag[1] = -952.38095238095241*dwdx2 - 952.38095238095241*dwdx3;
    JDiag[2] = -952.38095238095241*dwdx5 - 952.38095238095241*dwdx6;
    JDiag[3] = -952.38095238095241*dwdx7 - 952.38095238095241*dwdx8;
    JDiag[4] = -2857.1428571428573*dwdx9;
    JDiag[5] = -952.38095238095241*dwdx10 - 952.38095238095241*dwdx11;
    JDiag[6] = -2857.1428571428573*dwdx12;
    JDiag[7] = -952.38095238095241*dwdx13;
    JDiag[8] = -2857.1428571428573*dwdx14;
    JDiag[9] = -952.38095238095241*dwdx16;
    JDiag[10] = -952.38095238095241*dwdx17 - 952.38095238095241*dwdx18;
    JDiag[11] = -952.38095238095241*dwdx19 - 952.38095238095241*dwdx20 - 952.38095238095241*dwdx21;
    JDiag[12] = -952.38095238095241*dwdx22;
    JDiag[13] = -952.38095238095241*dwdx23 - 952.38095238095241*dwdx24;
    JDiag[14] = -952.38095238095241*dwdx25 - 952.38095238095241*dwdx26 - 952.38095238095241*dwdx27;
    JDiag[15] = -1.0*dwdx28;
}