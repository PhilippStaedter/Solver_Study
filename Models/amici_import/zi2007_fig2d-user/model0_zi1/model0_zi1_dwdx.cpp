#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_zi1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 0.0010499999999999999*kr_Cave;
    dwdx[1] = 0.0010499999999999999*Klid*Smads_Complex_n;
    dwdx[2] = 0.0010499999999999999*kr_EE;
    dwdx[3] = 0.0010499999999999999*Kcd;
    dwdx[4] = 0.0010499999999999999*Smad2c*Smad4c*k_Smads_Complex_c;
    dwdx[5] = 0.0010499999999999999*ki_Cave;
    dwdx[6] = 0.0010499999999999999*ki_EE;
    dwdx[7] = 0.0010499999999999999*Kimp_Smad2c;
    dwdx[8] = 0.0010499999999999999*LRC_EE*Smad4c*k_Smads_Complex_c;
    dwdx[9] = 0.00035*Kexp_Smad2n;
    dwdx[10] = 0.0010499999999999999*LRC_EE*Smad2c*k_Smads_Complex_c;
    dwdx[11] = 0.0010499999999999999*Kimp_Smad4c;
    dwdx[12] = 0.00035*Kexp_Smad4n;
    dwdx[13] = 0.0010499999999999999*Kimp_Smads_Complex_c;
    dwdx[14] = 0.00035*Kdiss_Smads_Complex_n;
    dwdx[15] = 0.0010499999999999999*Klid*LRC_Cave;
    dwdx[16] = 0.0010499999999999999*kr_Cave;
    dwdx[17] = 0.0010499999999999999*Kdeg_T1R_EE;
    dwdx[18] = 0.0010499999999999999*kr_EE;
    dwdx[19] = 0.0010499999999999999*T2R_Surf*TGF_beta*k_LRC;
    dwdx[20] = 0.0010499999999999999*ki_Cave;
    dwdx[21] = 0.0010499999999999999*ki_EE;
    dwdx[22] = 0.0010499999999999999*kr_Cave;
    dwdx[23] = 0.0010499999999999999*kr_EE;
    dwdx[24] = 0.0010499999999999999*Kdeg_T2R_EE;
    dwdx[25] = 0.0010499999999999999*ki_Cave;
    dwdx[26] = 0.0010499999999999999*ki_EE;
    dwdx[27] = 0.0010499999999999999*T1R_Surf*TGF_beta*k_LRC;
    dwdx[28] = 0.0010499999999999999*T1R_Surf*T2R_Surf*k_LRC;
}