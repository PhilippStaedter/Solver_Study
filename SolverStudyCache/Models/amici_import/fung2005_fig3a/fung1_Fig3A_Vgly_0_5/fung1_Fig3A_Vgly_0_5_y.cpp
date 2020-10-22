#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_fung1_Fig3A_Vgly_0_5(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = AcCoA;
    y[1] = AcP;
    y[2] = Acs;
    y[3] = HOAc;
    y[4] = HOAc_E;
    y[5] = LacI;
    y[6] = OAc;
    y[7] = Pta;
}