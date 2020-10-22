#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_bucher1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = ASL_b;
    y[1] = ASL_c;
    y[2] = ASL_m;
    y[3] = ASLoOH_b;
    y[4] = ASLoOH_c;
    y[5] = ASLoOH_m;
    y[6] = ASLpOH_b;
    y[7] = ASLpOH_c;
    y[8] = ASLpOH_m;
    y[9] = AS_b;
    y[10] = AS_c;
    y[11] = AS_m;
    y[12] = ASoOH_b;
    y[13] = ASoOH_c;
    y[14] = ASoOH_m;
    y[15] = ASpOH_b;
    y[16] = ASpOH_c;
    y[17] = ASpOH_m;
}