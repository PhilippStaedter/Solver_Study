#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_marhl(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*pow(Ca_cyt, 2)*v1_Kch*(CaER - Ca_cyt)/(pow(Ca_cyt, 2) + pow(v1_K1, 2));
    w[1] = 1.0*CaPr*v11_Kminus;
    w[2] = 1.0*Ca_cyt*Pr*v12_Kplus;
    w[3] = 1.0*v3_Kleak*(CaER - Ca_cyt);
    w[4] = 1.0*Ca_cyt*v5_Kpump;
    w[5] = 1.0*CaM*(pow(Ca_cyt, 2)*v7_Kout/(pow(Ca_cyt, 2) + pow(v7_K3, 2)) + v7_Km);
    w[6] = 1.0*pow(Ca_cyt, 8)*v9_Kin/(pow(Ca_cyt, 8) + pow(v9_K2, 8));
}