#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Leber2015(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Cdiff;
    y[1] = Commensal_Beneficial;
    y[2] = Commensal_Dead;
    y[3] = tDC_LP;
    y[4] = tDC_MLN;
    y[5] = Commensal_Harmful;
    y[6] = N_Lum;
    y[7] = E;
    y[8] = E_d;
    y[9] = iDC_E;
    y[10] = E_i;
    y[11] = M_LP;
    y[12] = eDC_LP;
    y[13] = M0;
    y[14] = N_LP;
    y[15] = Th17_LP;
    y[16] = Th1_LP;
    y[17] = iTreg_LP;
    y[18] = eDC_MLN;
    y[19] = iTreg_MLN;
    y[20] = nT;
    y[21] = Th17_MLN;
    y[22] = Th1_MLN;
}