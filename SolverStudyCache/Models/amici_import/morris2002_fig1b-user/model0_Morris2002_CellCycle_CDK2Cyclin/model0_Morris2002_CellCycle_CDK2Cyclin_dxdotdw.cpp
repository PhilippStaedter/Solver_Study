#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dxdotdw_model0_Morris2002_CellCycle_CDK2Cyclin(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdotdw[0] = -1000000000000.0;
    dxdotdw[1] = 1000000000000.0;
    dxdotdw[2] = 1000000000000.0;
    dxdotdw[3] = -1000000000000.0;
    dxdotdw[4] = -1000000000000.0;
}