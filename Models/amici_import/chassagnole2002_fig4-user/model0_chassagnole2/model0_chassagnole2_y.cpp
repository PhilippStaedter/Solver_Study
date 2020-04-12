#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_chassagnole2(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = cdhap;
    y[1] = ce4p;
    y[2] = cf6p;
    y[3] = cfdp;
    y[4] = cg1p;
    y[5] = cg6p;
    y[6] = cgap;
    y[7] = cglcex;
    y[8] = cpep;
    y[9] = cpg;
    y[10] = cpg2;
    y[11] = cpg3;
    y[12] = cpgp;
    y[13] = cpyr;
    y[14] = crib5p;
    y[15] = cribu5p;
    y[16] = csed7p;
    y[17] = cxyl5p;
}