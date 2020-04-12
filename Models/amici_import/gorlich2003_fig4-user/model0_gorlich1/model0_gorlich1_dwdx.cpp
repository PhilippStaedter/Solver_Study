#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_gorlich1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -1.2000000000000001e-11*GDP_dissociation_r7*RCC1_Ran;
    dwdx[1] = 1.2000000000000001e-11*GTP_binding_r3*RCC1_Ran;
    dwdx[2] = 1.2000000000000001e-11*RCC1_binding_r1*RanGDP_nuc;
    dwdx[3] = -1.2000000000000001e-11*RanGTP_nuc*RanGTP_release_r5;
    dwdx[4] = -1.2000000000000001e-11*GDP*GDP_dissociation_r7;
    dwdx[5] = 1.2000000000000001e-11*GTP*GTP_binding_r3;
    dwdx[6] = 1.2000000000000001e-11*GDP_dissociation_r2;
    dwdx[7] = -1.2000000000000001e-11*RCC1_binding_r8;
    dwdx[8] = -1.2000000000000001e-11*GTP_binding_r6;
    dwdx[9] = 1.2000000000000001e-11*RanGTP_release_r4;
    dwdx[10] = 1.7999999999999999e-11*RanGTP_binding_kon*RanGTP_cy;
    dwdx[11] = 1.7999999999999999e-11*RanBP1_RanGDP_kcat*RanGTP_RanBP1/(RanBP1_RanGDP_Km + RanGTP_RanBP1);
    dwdx[12] = 1.7999999999999999e-11*RanGAP_RanGDP_kcat_GAP*RanGTP_cy/(RanGAP_RanGDP_Km_GAP + RanGTP_cy);
    dwdx[13] = -1.2000000000000001e-11*Nucleoplasmic_transfer_kpermRanGDP;
    dwdx[14] = 1.2000000000000001e-11*Nucleoplasmic_transfer_kpermRanGDP;
    dwdx[15] = 1.2000000000000001e-11*RCC1*RCC1_binding_r1;
    dwdx[16] = -1.7999999999999999e-11*RanBP1_RanGDP_kcat*RanGAP*RanGTP_RanBP1/pow(RanBP1_RanGDP_Km + RanGTP_RanBP1, 2) + 1.7999999999999999e-11*RanBP1_RanGDP_kcat*RanGAP/(RanBP1_RanGDP_Km + RanGTP_RanBP1);
    dwdx[17] = -1.7999999999999999e-11*RanGTP_binding_koff;
    dwdx[18] = -1.2000000000000001e-11*Cytoplasmic_transfer_kpermRanGTP;
    dwdx[19] = -1.7999999999999999e-11*RanGAP*RanGAP_RanGDP_kcat_GAP*RanGTP_cy/pow(RanGAP_RanGDP_Km_GAP + RanGTP_cy, 2) + 1.7999999999999999e-11*RanGAP*RanGAP_RanGDP_kcat_GAP/(RanGAP_RanGDP_Km_GAP + RanGTP_cy);
    dwdx[20] = 1.7999999999999999e-11*RanBP1*RanGTP_binding_kon;
    dwdx[21] = 1.2000000000000001e-11*Cytoplasmic_transfer_kpermRanGTP;
    dwdx[22] = -1.2000000000000001e-11*RCC1*RanGTP_release_r5;
}