#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_orfao1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[2] = -1.0*PT;
            break;
        case 1:
            dwdp[2] = 1.0*PL*Va*Xa;
            break;
        case 2:
            dwdp[5] = 1.0*II*Xa/(II + km_2);
            break;
        case 3:
            dwdp[4] = 1.0*II*PT/(II + km_II);
            break;
        case 4:
            dwdp[3] = 1.0*IIa*V/(V + km_V);
            break;
        case 5:
            dwdp[0] = 1.0*RVV*X/(X + km_X);
            break;
        case 6:
            dwdp[7] = 1.0*IIa;
            break;
        case 7:
            dwdp[6] = 1.0*IIa;
            break;
        case 8:
            dwdp[1] = 1.0*Xa;
            break;
        case 9:
            dwdp[5] = -1.0*II*Xa*kcat_2/pow(II + km_2, 2);
            break;
        case 10:
            dwdp[4] = -1.0*II*PT*kcat_II/pow(II + km_II, 2);
            break;
        case 11:
            dwdp[3] = -1.0*IIa*V*kcat_V/pow(V + km_V, 2);
            break;
        case 12:
            dwdp[0] = -1.0*RVV*X*kcat_X/pow(X + km_X, 2);
            break;
    }
}