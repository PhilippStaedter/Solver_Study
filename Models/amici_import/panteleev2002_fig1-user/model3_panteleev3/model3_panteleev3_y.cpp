#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model3_panteleev3(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = TFPI;
    y[1] = VIIa_TF;
    y[2] = VIIa_TF_X;
    y[3] = VIIa_TF_Xa;
    y[4] = X;
    y[5] = Xa;
    y[6] = Xa_TFPI;
    y[7] = Xa_TFPI_VIIa_TF;
}