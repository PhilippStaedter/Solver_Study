#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_zi1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 0.0010499999999999999*Kdeg_T1R_EE*T1R_EE;
    w[1] = 0.0010499999999999999*v_T2R;
    w[2] = 0.0010499999999999999*T2R_Surf*ki_Cave;
    w[3] = 0.0010499999999999999*T2R_Cave*kr_Cave;
    w[4] = 0.0010499999999999999*T2R_Surf*ki_EE;
    w[5] = 0.0010499999999999999*T2R_EE*kr_EE;
    w[6] = 0.0010499999999999999*Kdeg_T2R_EE*T2R_EE;
    w[7] = 0.0010499999999999999*T1R_Surf*T2R_Surf*TGF_beta*k_LRC;
    w[8] = 0.0010499999999999999*LRC_Surf*ki_Cave;
    w[9] = 0.0010499999999999999*LRC_Cave*kr_Cave;
    w[10] = 0.0010499999999999999*Kimp_Smad2c*Smad2c;
    w[11] = 0.0010499999999999999*LRC_Surf*ki_EE;
    w[12] = 0.0010499999999999999*LRC_EE*kr_EE;
    w[13] = 0.0010499999999999999*Kcd*LRC_EE;
    w[14] = 0.0010499999999999999*LRC_EE*Smad2c*Smad4c*k_Smads_Complex_c;
    w[15] = 0.0010499999999999999*Kimp_Smads_Complex_c*Smads_Complex_c;
    w[16] = 0.00035*Kdiss_Smads_Complex_n*Smads_Complex_n;
    w[17] = 0.0010499999999999999*Klid*LRC_Cave*Smads_Complex_n;
    w[18] = 0.00035*Kexp_Smad2n*Smad2n;
    w[19] = 0.0010499999999999999*Kimp_Smad4c*Smad4c;
    w[20] = 0.00035*Kexp_Smad4n*Smad4n;
    w[21] = 0.0010499999999999999*v_T1R;
    w[22] = 0.0010499999999999999*T1R_Surf*ki_Cave;
    w[23] = 0.0010499999999999999*T1R_Cave*kr_Cave;
    w[24] = 0.0010499999999999999*T1R_Surf*ki_EE;
    w[25] = 0.0010499999999999999*T1R_EE*kr_EE;
}