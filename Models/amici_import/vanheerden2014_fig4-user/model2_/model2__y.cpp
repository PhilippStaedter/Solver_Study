#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model2_(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = EtOH;
    y[1] = Glycerol;
    y[2] = PiVac;
    y[3] = atp;
    y[4] = fbp;
    y[5] = g6p;
    y[6] = phos;
}