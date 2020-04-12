#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model2_(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[24] = -2.0*dwdx0 + 4.0*dwdx1 - 1.0*dwdx2;
    J[25] = 4.0*dwdx3;
    J[27] = 4.0*dwdx5;
    J[31] = 1.0*dwdx0 - 1.0*dwdx1;
    J[32] = -1.0*dwdx3 - 1.0*dwdx4;
    J[34] = -1.0*dwdx5;
    J[45] = -2.0*dwdx1 + 1.0*dwdx2;
    J[46] = -2.0*dwdx3 + 2.0*dwdx4;
    J[48] = -2.0*dwdx5 + 1.0*dwdx6;
}