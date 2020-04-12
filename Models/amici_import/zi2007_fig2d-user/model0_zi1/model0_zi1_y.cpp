#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_zi1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = LRC_Cave;
    y[1] = LRC_EE;
    y[2] = LRC_Surf;
    y[3] = Smad2c;
    y[4] = Smad2n;
    y[5] = Smad4c;
    y[6] = Smad4n;
    y[7] = Smads_Complex_c;
    y[8] = Smads_Complex_n;
    y[9] = T1R_Cave;
    y[10] = T1R_EE;
    y[11] = T1R_Surf;
    y[12] = T2R_Cave;
    y[13] = T2R_EE;
    y[14] = T2R_Surf;
    y[15] = TGF_beta;
}