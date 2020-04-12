#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Ueda2001(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = EmptySet;
    y[1] = CCc;
    y[2] = CCn;
    y[3] = Clkc;
    y[4] = Clkm;
    y[5] = Perc;
    y[6] = Perm;
    y[7] = PTc;
    y[8] = PTn;
    y[9] = Timc;
    y[10] = Timm;
    y[11] = species_0000012;
    y[12] = species_0000013;
}