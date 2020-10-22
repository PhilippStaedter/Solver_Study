#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_brands1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = AMP;
    y[1] = Acetic_acid;
    y[2] = Amadori;
    y[3] = C5;
    y[4] = Cn;
    y[5] = Formic_acid;
    y[6] = Fru;
    y[7] = Glu;
    y[8] = Melanoidin;
    y[9] = Triose;
    y[10] = lys_R;
}