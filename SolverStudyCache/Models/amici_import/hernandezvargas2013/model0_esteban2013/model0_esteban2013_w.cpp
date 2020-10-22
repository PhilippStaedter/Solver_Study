#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_esteban2013(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = sT;
    w[1] = MStar*deltaMStar;
    w[2] = TStar*pT;
    w[3] = MStar*pM;
    w[4] = V*deltaV;
    w[5] = T*V*rhoT/(CT + V);
    w[6] = T*V*kT;
    w[7] = T*deltaT;
    w[8] = TStar*deltaTStar;
    w[9] = sM;
    w[10] = M*V*rhoM/(CM + V);
    w[11] = M*V*kM;
    w[12] = M*deltaM;
}