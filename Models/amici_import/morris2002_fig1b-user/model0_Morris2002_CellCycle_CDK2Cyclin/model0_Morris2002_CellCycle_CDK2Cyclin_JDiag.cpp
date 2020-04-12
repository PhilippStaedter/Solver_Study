#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model0_Morris2002_CellCycle_CDK2Cyclin(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1000000000000.0*dwdx0 + 1000000000000.0*dwdx1;
    JDiag[1] = 1000000000000.0*dwdx2;
    JDiag[2] = -1000000000000.0*dwdx3;
    JDiag[3] = -1000000000000.0*dwdx4;
}