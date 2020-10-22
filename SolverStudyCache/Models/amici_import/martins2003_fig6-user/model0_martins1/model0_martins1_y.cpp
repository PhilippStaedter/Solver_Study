#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_martins1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = AA;
    y[1] = Cn;
    y[2] = DFG;
    y[3] = E1;
    y[4] = E2;
    y[5] = FA;
    y[6] = Fru;
    y[7] = Glu;
    y[8] = Gly;
    y[9] = MG;
    y[10] = Man;
    y[11] = Mel;
    y[12] = _1DG;
    y[13] = _3DG;
}