#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_orfao1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*RVV*X*kcat_X/(X + km_X);
    w[1] = 1.0*Xa*ki_Xa;
    w[2] = 1.0*PL*Va*Xa*k_PT - 1.0*PT*k_PL;
    w[3] = 1.0*IIa*V*kcat_V/(V + km_V);
    w[4] = 1.0*II*PT*kcat_II/(II + km_II);
    w[5] = 1.0*II*Xa*kcat_2/(II + km_2);
    w[6] = 1.0*IIa*ki_IIaAlpha2M;
    w[7] = 1.0*IIa*ki_IIaATIII;
}