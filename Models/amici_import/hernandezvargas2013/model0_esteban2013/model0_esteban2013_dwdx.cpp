#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_esteban2013(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = V*rhoM/(CM + V);
    dwdx[1] = V*kM;
    dwdx[2] = deltaM;
    dwdx[3] = deltaMStar;
    dwdx[4] = pM;
    dwdx[5] = V*rhoT/(CT + V);
    dwdx[6] = V*kT;
    dwdx[7] = deltaT;
    dwdx[8] = pT;
    dwdx[9] = deltaTStar;
    dwdx[10] = deltaV;
    dwdx[11] = -T*V*rhoT/pow(CT + V, 2) + T*rhoT/(CT + V);
    dwdx[12] = T*kT;
    dwdx[13] = -M*V*rhoM/pow(CM + V, 2) + M*rhoM/(CM + V);
    dwdx[14] = M*kM;
}