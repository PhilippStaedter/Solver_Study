#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_Leber2015(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 0.070000000000000007*Treg_Degradation_k1*iTreg_LP;
    w[1] = 1.0*eDC_Degradation_k1*eDC_MLN;
    w[2] = 0.070000000000000007*Th17_Degradation_k1*Th17_LP;
    w[3] = 0.070000000000000007*Th1_Degradation_k1*Th1_LP;
    w[4] = 1.0*Commensal_Beneficial*N_Degradation_K*N_Lum;
    w[5] = 4.0*E*E_Damage_v*(E_Damage_k1*N_Lum + E_Damage_k2*Th17_LP + E_Damage_k3*M_LP);
    w[6] = eDC_LP*eDC_Migration_k1;
    w[7] = Cdiff*eDC_Production_k;
    w[8] = 1.0*Cdiff*Cdiff_Death_K*(Cdiff_Death_m2*N_Lum - Cdiff_Death_m3*Commensal_Harmful + M_LP);
    w[9] = N_Activation_Migration_v*N_LP*(Cdiff*(E_d*N_Activation_Migration_k1 + N_Activation_Migration_k2*Th17_LP) - N_Activation_Migration_k3*iTreg_LP);
    w[10] = 1.0*Cdiff*Cdiff_Growth_K*Commensal_Harmful;
    w[11] = Treg_Migration_k1*iTreg_MLN;
    w[12] = Th1_MLN*Th1_Migration_k1;
    w[13] = -0.070000000000000007*Cdiff*Th17_Plasticity_k2*iTreg_LP + 0.070000000000000007*Th17_LP*Th17_Plasticity_k1;
    w[14] = Th17_MLN*Th17_Migration_k1;
    w[15] = 4.0*Cdiff*E*E_Inflame_K;
    w[16] = 4.0*E_i*E_i_Damage_v*(E_i_Damage_k1*N_Lum + E_i_Damage_k2*Th17_LP + E_i_Damage_k3*M_LP);
    w[17] = M0*M_Activation_K*(Cdiff + M_Activation_e1*Th17_LP - M_Activation_e2*iTreg_LP);
    w[18] = 4.0*M_Death_k1*M_LP;
    w[19] = 1.0*Commensal_Beneficial*Commensal_Regrowth_k1*E_i*N_Lum - 1.0*Commensal_Dead*Commensal_Regrowth_k2;
    w[20] = 4.0*E_Heal_k1*E_d;
    w[21] = Cdiff*tDC_Production_K*(Commensal_Beneficial*tDC_Production_k1/Commensal_Dead + E*tDC_Production_k2/(E_i + 100));
    w[22] = 1.0*tDC_LP*tDC_Migration_k1;
    w[23] = 1.0*tDC_Degradation_k*tDC_MLN;
    w[24] = 1.0*Th17_Differentiation_k1*eDC_MLN;
    w[25] = 1.0*Commensal_Dead*Th1_Differentiation_K*eDC_MLN/(Commensal_Beneficial*Th1_Differentiation_k1 + E*Th1_Differentiation_k2);
    w[26] = Treg_Differentiation_k1*tDC_MLN;
    w[27] = 1.0*Commensal_Harmful*Commensal_Harmful_Death_K*(Commensal_Harmful_Death_A1*N_LP + Commensal_Harmful_Death_A2*E_i);
    w[28] = 1.0*Commensal_Dead*Commensal_Death_k1;
    w[29] = 4.0*E_i*E_i_Natural_Death_k1;
}