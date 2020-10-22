#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Leber2015(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 0.070000000000000007*iTreg_LP;
            break;
        case 1:
            dwdp[1] = 1.0*eDC_MLN;
            break;
        case 2:
            dwdp[2] = 0.070000000000000007*Th17_LP;
            break;
        case 3:
            dwdp[3] = 0.070000000000000007*Th1_LP;
            break;
        case 4:
            dwdp[4] = 1.0*Commensal_Beneficial*N_Lum;
            break;
        case 5:
            dwdp[5] = 4.0*E*E_Damage_v*M_LP;
            break;
        case 6:
            dwdp[5] = 4.0*E*E_Damage_v*Th17_LP;
            break;
        case 7:
            dwdp[5] = 4.0*E*E_Damage_v*N_Lum;
            break;
        case 8:
            dwdp[5] = 4.0*E*(E_Damage_k1*N_Lum + E_Damage_k2*Th17_LP + E_Damage_k3*M_LP);
            break;
        case 9:
            dwdp[6] = eDC_LP;
            break;
        case 10:
            dwdp[7] = Cdiff;
            break;
        case 11:
            dwdp[8] = -1.0*Cdiff*Cdiff_Death_K*Commensal_Harmful;
            break;
        case 12:
            dwdp[8] = 1.0*Cdiff*Cdiff_Death_K*N_Lum;
            break;
        case 13:
            dwdp[8] = 1.0*Cdiff*(Cdiff_Death_m2*N_Lum - Cdiff_Death_m3*Commensal_Harmful + M_LP);
            break;
        case 14:
            dwdp[9] = -N_Activation_Migration_v*N_LP*iTreg_LP;
            break;
        case 15:
            dwdp[9] = Cdiff*N_Activation_Migration_v*N_LP*Th17_LP;
            break;
        case 16:
            dwdp[9] = Cdiff*E_d*N_Activation_Migration_v*N_LP;
            break;
        case 17:
            dwdp[9] = N_LP*(Cdiff*(E_d*N_Activation_Migration_k1 + N_Activation_Migration_k2*Th17_LP) - N_Activation_Migration_k3*iTreg_LP);
            break;
        case 18:
            dwdp[10] = 1.0*Cdiff*Commensal_Harmful;
            break;
        case 19:
            dwdp[11] = iTreg_MLN;
            break;
        case 20:
            dwdp[12] = Th1_MLN;
            break;
        case 21:
            dwdp[13] = -0.070000000000000007*Cdiff*iTreg_LP;
            break;
        case 22:
            dwdp[13] = 0.070000000000000007*Th17_LP;
            break;
        case 23:
            dwdp[14] = Th17_MLN;
            break;
        case 24:
            dwdp[15] = 4.0*Cdiff*E;
            break;
        case 25:
            dwdp[16] = 4.0*E_i*E_i_Damage_v*M_LP;
            break;
        case 26:
            dwdp[16] = 4.0*E_i*E_i_Damage_v*Th17_LP;
            break;
        case 27:
            dwdp[16] = 4.0*E_i*E_i_Damage_v*N_Lum;
            break;
        case 28:
            dwdp[16] = 4.0*E_i*(E_i_Damage_k1*N_Lum + E_i_Damage_k2*Th17_LP + E_i_Damage_k3*M_LP);
            break;
        case 29:
            dwdp[17] = -M0*M_Activation_K*iTreg_LP;
            break;
        case 30:
            dwdp[17] = M0*M_Activation_K*Th17_LP;
            break;
        case 31:
            dwdp[17] = M0*(Cdiff + M_Activation_e1*Th17_LP - M_Activation_e2*iTreg_LP);
            break;
        case 32:
            dwdp[18] = 4.0*M_LP;
            break;
        case 33:
            dwdp[19] = -1.0*Commensal_Dead;
            break;
        case 34:
            dwdp[19] = 1.0*Commensal_Beneficial*E_i*N_Lum;
            break;
        case 35:
            dwdp[20] = 4.0*E_d;
            break;
        case 36:
            dwdp[21] = Cdiff*E*tDC_Production_K/(E_i + 100);
            break;
        case 37:
            dwdp[21] = Cdiff*Commensal_Beneficial*tDC_Production_K/Commensal_Dead;
            break;
        case 38:
            dwdp[21] = Cdiff*(Commensal_Beneficial*tDC_Production_k1/Commensal_Dead + E*tDC_Production_k2/(E_i + 100));
            break;
        case 39:
            dwdp[22] = 1.0*tDC_LP;
            break;
        case 40:
            dwdp[23] = 1.0*tDC_MLN;
            break;
        case 41:
            dwdp[24] = 1.0*eDC_MLN;
            break;
        case 42:
            dwdp[25] = -1.0*Commensal_Beneficial*Commensal_Dead*Th1_Differentiation_K*eDC_MLN/pow(Commensal_Beneficial*Th1_Differentiation_k1 + E*Th1_Differentiation_k2, 2);
            break;
        case 43:
            dwdp[25] = -1.0*Commensal_Dead*E*Th1_Differentiation_K*eDC_MLN/pow(Commensal_Beneficial*Th1_Differentiation_k1 + E*Th1_Differentiation_k2, 2);
            break;
        case 44:
            dwdp[25] = 1.0*Commensal_Dead*eDC_MLN/(Commensal_Beneficial*Th1_Differentiation_k1 + E*Th1_Differentiation_k2);
            break;
        case 45:
            dwdp[26] = tDC_MLN;
            break;
        case 46:
            dwdp[27] = 1.0*Commensal_Harmful*Commensal_Harmful_Death_K*E_i;
            break;
        case 47:
            dwdp[27] = 1.0*Commensal_Harmful*Commensal_Harmful_Death_K*N_LP;
            break;
        case 48:
            dwdp[27] = 1.0*Commensal_Harmful*(Commensal_Harmful_Death_A1*N_LP + Commensal_Harmful_Death_A2*E_i);
            break;
        case 49:
            dwdp[28] = 1.0*Commensal_Dead;
            break;
        case 50:
            dwdp[29] = 4.0*E_i;
            break;
    }
}