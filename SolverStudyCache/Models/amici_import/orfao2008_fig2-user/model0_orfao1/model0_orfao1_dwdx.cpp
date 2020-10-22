#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_orfao1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -1.0*II*PT*kcat_II/pow(II + km_II, 2) + 1.0*PT*kcat_II/(II + km_II);
    dwdx[1] = -1.0*II*Xa*kcat_2/pow(II + km_2, 2) + 1.0*Xa*kcat_2/(II + km_2);
    dwdx[2] = 1.0*V*kcat_V/(V + km_V);
    dwdx[3] = 1.0*ki_IIaAlpha2M;
    dwdx[4] = 1.0*ki_IIaATIII;
    dwdx[5] = 1.0*Va*Xa*k_PT;
    dwdx[6] = -1.0*k_PL;
    dwdx[7] = 1.0*II*kcat_II/(II + km_II);
    dwdx[8] = 1.0*X*kcat_X/(X + km_X);
    dwdx[9] = -1.0*IIa*V*kcat_V/pow(V + km_V, 2) + 1.0*IIa*kcat_V/(V + km_V);
    dwdx[10] = 1.0*PL*Xa*k_PT;
    dwdx[11] = -1.0*RVV*X*kcat_X/pow(X + km_X, 2) + 1.0*RVV*kcat_X/(X + km_X);
    dwdx[12] = 1.0*ki_Xa;
    dwdx[13] = 1.0*PL*Va*k_PT;
    dwdx[14] = 1.0*II*kcat_2/(II + km_2);
}