#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_marhl(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*pow(Ca_cyt, 2)*v1_Kch/(pow(Ca_cyt, 2) + pow(v1_K1, 2));
    dwdx[1] = 1.0*v3_Kleak;
    dwdx[2] = 1.0*pow(Ca_cyt, 2)*v7_Kout/(pow(Ca_cyt, 2) + pow(v7_K3, 2)) + 1.0*v7_Km;
    dwdx[3] = 1.0*v11_Kminus;
    dwdx[4] = -2.0*pow(Ca_cyt, 3)*v1_Kch*(CaER - Ca_cyt)/pow(pow(Ca_cyt, 2) + pow(v1_K1, 2), 2) - 1.0*pow(Ca_cyt, 2)*v1_Kch/(pow(Ca_cyt, 2) + pow(v1_K1, 2)) + 2.0*Ca_cyt*v1_Kch*(CaER - Ca_cyt)/(pow(Ca_cyt, 2) + pow(v1_K1, 2));
    dwdx[5] = 1.0*Pr*v12_Kplus;
    dwdx[6] = -1.0*v3_Kleak;
    dwdx[7] = 1.0*v5_Kpump;
    dwdx[8] = 1.0*CaM*(-2*pow(Ca_cyt, 3)*v7_Kout/pow(pow(Ca_cyt, 2) + pow(v7_K3, 2), 2) + 2*Ca_cyt*v7_Kout/(pow(Ca_cyt, 2) + pow(v7_K3, 2)));
    dwdx[9] = -8.0*pow(Ca_cyt, 15)*v9_Kin/pow(pow(Ca_cyt, 8) + pow(v9_K2, 8), 2) + 8.0*pow(Ca_cyt, 7)*v9_Kin/(pow(Ca_cyt, 8) + pow(v9_K2, 8));
    dwdx[10] = 1.0*Ca_cyt*v12_Kplus;
}