#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model2_(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[24] = 2.0*dwdx0 - 4.0*dwdx1 + 1.0*dwdx2;
    JB[25] = -1.0*dwdx0 + 1.0*dwdx1;
    JB[27] = 2.0*dwdx1 - 1.0*dwdx2;
    JB[31] = -4.0*dwdx3;
    JB[32] = 1.0*dwdx3 + 1.0*dwdx4;
    JB[34] = 2.0*dwdx3 - 2.0*dwdx4;
    JB[45] = -4.0*dwdx5;
    JB[46] = 1.0*dwdx5;
    JB[48] = 2.0*dwdx5 - 1.0*dwdx6;
}