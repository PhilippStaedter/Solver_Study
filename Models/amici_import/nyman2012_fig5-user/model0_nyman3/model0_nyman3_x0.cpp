#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_nyman3(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 8.9406759753263199;
    x0[1] = 9.4399819422554394;
    x0[2] = 0.56001805774457303;
    x0[3] = 4.8386389075851502e-6;
    x0[4] = 0.42407663182338401;
    x0[5] = 0.59688996214639001;
    x0[6] = 0.0383525925240207;
    x0[7] = 9.9963588640715102;
    x0[8] = 0.0036411359284838599;
}