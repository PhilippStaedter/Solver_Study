#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_ferreira1_Fig2F(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Amadori;
    y[1] = CML;
    y[2] = Glucose;
    y[3] = Glyoxal;
    y[4] = Lysine;
    y[5] = Schiff;
}