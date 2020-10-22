#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_leloup1_Fig2A(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = CC;
    y[1] = Cn;
    y[2] = Mp;
    y[3] = Mt;
    y[4] = P0;
    y[5] = P1;
    y[6] = P2;
    y[7] = T0;
    y[8] = T1;
    y[9] = T2;
}