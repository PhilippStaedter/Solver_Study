#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_Leber2015(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = eDC_Production_k;
    dwdx[1] = 1.0*Cdiff_Death_K*(Cdiff_Death_m2*N_Lum - Cdiff_Death_m3*Commensal_Harmful + M_LP);
    dwdx[2] = N_Activation_Migration_v*N_LP*(E_d*N_Activation_Migration_k1 + N_Activation_Migration_k2*Th17_LP);
    dwdx[3] = 1.0*Cdiff_Growth_K*Commensal_Harmful;
    dwdx[4] = -0.070000000000000007*Th17_Plasticity_k2*iTreg_LP;
    dwdx[5] = 4.0*E*E_Inflame_K;
    dwdx[6] = M0*M_Activation_K;
    dwdx[7] = tDC_Production_K*(Commensal_Beneficial*tDC_Production_k1/Commensal_Dead + E*tDC_Production_k2/(E_i + 100));
    dwdx[8] = 1.0*N_Degradation_K*N_Lum;
    dwdx[9] = 1.0*Commensal_Regrowth_k1*E_i*N_Lum;
    dwdx[10] = Cdiff*tDC_Production_K*tDC_Production_k1/Commensal_Dead;
    dwdx[11] = -1.0*Commensal_Dead*Th1_Differentiation_K*Th1_Differentiation_k1*eDC_MLN/pow(Commensal_Beneficial*Th1_Differentiation_k1 + E*Th1_Differentiation_k2, 2);
    dwdx[12] = -1.0*Commensal_Regrowth_k2;
    dwdx[13] = -Cdiff*Commensal_Beneficial*tDC_Production_K*tDC_Production_k1/pow(Commensal_Dead, 2);
    dwdx[14] = 1.0*Th1_Differentiation_K*eDC_MLN/(Commensal_Beneficial*Th1_Differentiation_k1 + E*Th1_Differentiation_k2);
    dwdx[15] = 1.0*Commensal_Death_k1;
    dwdx[16] = 1.0*tDC_Migration_k1;
    dwdx[17] = 1.0*tDC_Degradation_k;
    dwdx[18] = Treg_Differentiation_k1;
    dwdx[19] = -1.0*Cdiff*Cdiff_Death_K*Cdiff_Death_m3;
    dwdx[20] = 1.0*Cdiff*Cdiff_Growth_K;
    dwdx[21] = 1.0*Commensal_Harmful_Death_K*(Commensal_Harmful_Death_A1*N_LP + Commensal_Harmful_Death_A2*E_i);
    dwdx[22] = 1.0*Commensal_Beneficial*N_Degradation_K;
    dwdx[23] = 4.0*E*E_Damage_k1*E_Damage_v;
    dwdx[24] = 1.0*Cdiff*Cdiff_Death_K*Cdiff_Death_m2;
    dwdx[25] = 4.0*E_i*E_i_Damage_k1*E_i_Damage_v;
    dwdx[26] = 1.0*Commensal_Beneficial*Commensal_Regrowth_k1*E_i;
    dwdx[27] = 4.0*E_Damage_v*(E_Damage_k1*N_Lum + E_Damage_k2*Th17_LP + E_Damage_k3*M_LP);
    dwdx[28] = 4.0*Cdiff*E_Inflame_K;
    dwdx[29] = Cdiff*tDC_Production_K*tDC_Production_k2/(E_i + 100);
    dwdx[30] = -1.0*Commensal_Dead*Th1_Differentiation_K*Th1_Differentiation_k2*eDC_MLN/pow(Commensal_Beneficial*Th1_Differentiation_k1 + E*Th1_Differentiation_k2, 2);
    dwdx[31] = Cdiff*N_Activation_Migration_k1*N_Activation_Migration_v*N_LP;
    dwdx[32] = 4.0*E_Heal_k1;
    dwdx[33] = 4.0*E_i_Damage_v*(E_i_Damage_k1*N_Lum + E_i_Damage_k2*Th17_LP + E_i_Damage_k3*M_LP);
    dwdx[34] = 1.0*Commensal_Beneficial*Commensal_Regrowth_k1*N_Lum;
    dwdx[35] = -Cdiff*E*tDC_Production_K*tDC_Production_k2/pow(E_i + 100, 2);
    dwdx[36] = 1.0*Commensal_Harmful*Commensal_Harmful_Death_A2*Commensal_Harmful_Death_K;
    dwdx[37] = 4.0*E_i_Natural_Death_k1;
    dwdx[38] = 4.0*E*E_Damage_k3*E_Damage_v;
    dwdx[39] = 1.0*Cdiff*Cdiff_Death_K;
    dwdx[40] = 4.0*E_i*E_i_Damage_k3*E_i_Damage_v;
    dwdx[41] = 4.0*M_Death_k1;
    dwdx[42] = eDC_Migration_k1;
    dwdx[43] = M_Activation_K*(Cdiff + M_Activation_e1*Th17_LP - M_Activation_e2*iTreg_LP);
    dwdx[44] = N_Activation_Migration_v*(Cdiff*(E_d*N_Activation_Migration_k1 + N_Activation_Migration_k2*Th17_LP) - N_Activation_Migration_k3*iTreg_LP);
    dwdx[45] = 1.0*Commensal_Harmful*Commensal_Harmful_Death_A1*Commensal_Harmful_Death_K;
    dwdx[46] = 0.070000000000000007*Th17_Degradation_k1;
    dwdx[47] = 4.0*E*E_Damage_k2*E_Damage_v;
    dwdx[48] = Cdiff*N_Activation_Migration_k2*N_Activation_Migration_v*N_LP;
    dwdx[49] = 0.070000000000000007*Th17_Plasticity_k1;
    dwdx[50] = 4.0*E_i*E_i_Damage_k2*E_i_Damage_v;
    dwdx[51] = M0*M_Activation_K*M_Activation_e1;
    dwdx[52] = 0.070000000000000007*Th1_Degradation_k1;
    dwdx[53] = 0.070000000000000007*Treg_Degradation_k1;
    dwdx[54] = -N_Activation_Migration_k3*N_Activation_Migration_v*N_LP;
    dwdx[55] = -0.070000000000000007*Cdiff*Th17_Plasticity_k2;
    dwdx[56] = -M0*M_Activation_K*M_Activation_e2;
    dwdx[57] = 1.0*eDC_Degradation_k1;
    dwdx[58] = 1.0*Th17_Differentiation_k1;
    dwdx[59] = 1.0*Commensal_Dead*Th1_Differentiation_K/(Commensal_Beneficial*Th1_Differentiation_k1 + E*Th1_Differentiation_k2);
    dwdx[60] = Treg_Migration_k1;
    dwdx[61] = Th17_Migration_k1;
    dwdx[62] = Th1_Migration_k1;
}