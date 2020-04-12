#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model2_saeidi1(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0;
    JDiag[1] = -1.0*dwdx1;
    JDiag[2] = -1.0*dwdx2 + 1.0*dwdx3;
    JDiag[3] = -1.0*dwdx4 - 1.0*dwdx5;
    JDiag[4] = 1.0*dwdx6 - 1.0*dwdx7;
    JDiag[6] = -1.0*dwdx8 - 1.0*dwdx9;
    JDiag[7] = 1.0*dwdx10 - 1.0*dwdx11;
    JDiag[8] = 1.0*dwdx12;
}