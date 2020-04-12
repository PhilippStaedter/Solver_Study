#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model2_perelson3(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = N0*Tstarstar*muB;
    w[1] = T*V*k1;
    w[2] = Tstar*k2;
    w[3] = V*muVV;
    w[4] = Tstarstar*muB;
    w[5] = T*muT;
    w[6] = Tstar*muT;
    w[7] = T*r*(1 - Ttot/Tmax);
    w[8] = s*theta/(V + theta);
}