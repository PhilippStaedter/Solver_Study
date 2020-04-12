#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_miao2008(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = T*rho;
    w[1] = T*Tm*kM;
    w[2] = T*Tw*kW;
    w[3] = T*Tmw*kR;
    w[4] = Tm*rhoM;
    w[5] = Tm*Tw*qM;
    w[6] = Tw*rhoW;
    w[7] = Tm*Tw*qW;
    w[8] = Tmw*rhoMW;
}