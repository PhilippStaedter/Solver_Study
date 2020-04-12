#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_miao2008(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = rho;
    dwdx[1] = Tm*kM;
    dwdx[2] = Tw*kW;
    dwdx[3] = Tmw*kR;
    dwdx[4] = T*kM;
    dwdx[5] = rhoM;
    dwdx[6] = Tw*qM;
    dwdx[7] = Tw*qW;
    dwdx[8] = T*kR;
    dwdx[9] = rhoMW;
    dwdx[10] = T*kW;
    dwdx[11] = Tm*qM;
    dwdx[12] = rhoW;
    dwdx[13] = Tm*qW;
}