#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_kowald1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 0.01*k10*species_0000001;
    w[1] = 1.0*k2*species_0000001*species_0000002;
    w[2] = 1.0*k1;
    w[3] = 1.0*k3*species_0000001*(-species_0000002 + species_0000016);
    w[4] = 1.0*k13a*(-species_0000002 + species_0000016);
    w[5] = 1.0*k13b*species_0000002;
    w[6] = 1.0*k17*species_0000011;
    w[7] = 1.0*k18*species_0000007;
    w[8] = 1.0*k19*pow(species_0000007, 2);
    w[9] = 1.0*k4*species_0000001*species_0000007;
    w[10] = 1.0*k5*species_0000001*species_0000006;
    w[11] = 1.0*k6*species_0000002*species_0000006;
    w[12] = 1.0*k7*species_0000006*species_0000017;
    w[13] = 1.0*k9*species_0000008;
    w[14] = 0.01*k10*species_0000001;
    w[15] = 1.0*k11*species_0000008;
    w[16] = 1.0*k12*species_0000009;
}