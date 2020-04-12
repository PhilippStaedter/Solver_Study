#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_ferreira1_Fig2F(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = -1.0*dwdx0 - 1.0*dwdx1;
    JSparse[1] = 1.0*dwdx1;
    JSparse[2] = 1.0*dwdx0;
    JSparse[3] = -1.0*dwdx2 - 1.0*dwdx3;
    JSparse[4] = 1.0*dwdx3;
    JSparse[5] = -1.0*dwdx2;
    JSparse[6] = 1.0*dwdx2;
    JSparse[7] = 1.0*dwdx4;
    JSparse[8] = -1.0*dwdx4 - 1.0*dwdx5;
    JSparse[9] = -1.0*dwdx4;
    JSparse[10] = 1.0*dwdx7;
    JSparse[11] = -1.0*dwdx6;
    JSparse[12] = -1.0*dwdx7;
    JSparse[13] = -1.0*dwdx6 - 1.0*dwdx7;
    JSparse[14] = 1.0*dwdx6;
    JSparse[15] = 1.0*dwdx9;
    JSparse[16] = 1.0*dwdx10;
    JSparse[17] = 1.0*dwdx8;
    JSparse[18] = 1.0*dwdx12;
    JSparse[19] = 1.0*dwdx11 + 1.0*dwdx8;
    JSparse[20] = -1.0*dwdx10 - 1.0*dwdx13 - 1.0*dwdx8 - 1.0*dwdx9;
}