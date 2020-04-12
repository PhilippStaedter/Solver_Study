#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_brands1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = Glu*_J1_K1;
    w[1] = Fru*_J10_K10*lys_R;
    w[2] = AMP*_J11_K11;
    w[3] = Fru*_J2_K2;
    w[4] = Glu*_J3_K3;
    w[5] = Fru*_J4_K4;
    w[6] = Fru*_J5_K5;
    w[7] = Triose*_J6_K6;
    w[8] = Glu*_J7_K7*lys_R;
    w[9] = Amadori*_J8_K8;
    w[10] = Amadori*_J9_K9;
}