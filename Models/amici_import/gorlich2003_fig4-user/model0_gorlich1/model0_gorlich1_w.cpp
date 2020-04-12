#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_gorlich1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.2000000000000001e-11*Cytoplasmic_transfer_kpermRanGTP*(-RanGTP_cy + RanGTP_nuc);
    w[1] = -1.2000000000000001e-11*GDP*GDP_dissociation_r7*RCC1_Ran + 1.2000000000000001e-11*GDP_dissociation_r2*RCC1_RanGDP;
    w[2] = 1.2000000000000001e-11*GTP*GTP_binding_r3*RCC1_Ran - 1.2000000000000001e-11*GTP_binding_r6*RCC1_RanGTP;
    w[3] = 1.2000000000000001e-11*Nucleoplasmic_transfer_kpermRanGDP*(-RanGDP_cy + RanGDP_nuc);
    w[4] = 1.2000000000000001e-11*RCC1*RCC1_binding_r1*RanGDP_nuc - 1.2000000000000001e-11*RCC1_RanGDP*RCC1_binding_r8;
    w[5] = 1.7999999999999999e-11*RanBP1_RanGDP_kcat*RanGAP*RanGTP_RanBP1/(RanBP1_RanGDP_Km + RanGTP_RanBP1);
    w[6] = 1.7999999999999999e-11*RanGAP*RanGAP_RanGDP_kcat_GAP*RanGTP_cy/(RanGAP_RanGDP_Km_GAP + RanGTP_cy);
    w[7] = 1.7999999999999999e-11*RanBP1*RanGTP_binding_kon*RanGTP_cy - 1.7999999999999999e-11*RanGTP_RanBP1*RanGTP_binding_koff;
    w[8] = -1.2000000000000001e-11*RCC1*RanGTP_nuc*RanGTP_release_r5 + 1.2000000000000001e-11*RCC1_RanGTP*RanGTP_release_r4;
}