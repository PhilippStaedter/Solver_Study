#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_neves1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -(2.7593514420370556e-6*A1_Kf*c2_R2C2_cyto*cAMP_cyto - 0.0016611295681063099*A1_Kr*c3_R2C2_cyto)/pow(KMOLE, 2);
            dwdp[1] = -(2.7593514420370556e-6*A2_Kf*c3_R2C2_cyto*cAMP_cyto - 0.0016611295681063099*A2_Kr*PKA_cyto)/pow(KMOLE, 2);
            dwdp[5] = -(2.7593514420370556e-6*B1_Kf*R2C2_cyto*cAMP_cyto - 0.0016611295681063099*B1_Kr*c_R2C2_cyto)/pow(KMOLE, 2);
            dwdp[6] = -(2.7593514420370556e-6*B2_Kf*cAMP_cyto*c_R2C2_cyto - 0.0016611295681063099*B2_Kr*c2_R2C2_cyto)/pow(KMOLE, 2);
            dwdp[9] = -(0.0016611295681063099*GTPase_Kf_GTPase*G_a_s_cyto - 0.0016611295681063099*GTPase_Kr_GTPase*G_GDP_cyto)/pow(KMOLE, 2);
            dwdp[12] = -4.1390271630555837e-7*MAPK_cyto*MEK_active_cyto/(pow(KMOLE, 2)*(0.0016611295681063099*MAPK_cyto + MEK_activates_MAPK_Km));
            dwdp[13] = -2.2074811536296445e-5*PDE4_cyto*cAMP_cyto/(pow(KMOLE, 2)*(PDE4_Km_PDE4 + 0.0016611295681063099*cAMP_cyto));
            dwdp[14] = -2.759351442037056e-5*PDE4_cyto*PKA_cyto/(pow(KMOLE, 2)*(0.0016611295681063099*PDE4_cyto + PKA_P_PDE_Km));
            dwdp[15] = -2.7593514420370556e-6*PKA_cyto*PTP_cyto*kcat_PKA_P_PTP/(pow(KMOLE, 2)*(PKA_P_PTP_Km + 0.0016611295681063099*PTP_cyto));
            dwdp[16] = -2.7593514420370556e-6*B_Raf_cyto*PKA_cyto*kcat_PKA_activates_Raf/(pow(KMOLE, 2)*(0.0016611295681063099*B_Raf_cyto + PKA_activates_Raf_Km));
            dwdp[17] = -2.7593514420370556e-6*MAPK_active_cyto*PP2A_cyto*kcat_PPase_MAPK/(pow(KMOLE, 2)*(0.0016611295681063099*MAPK_active_cyto + PPase_MAPK_Km));
            dwdp[18] = -2.7593514420370556e-6*B_Raf_active_cyto*PP2A_cyto*kcat_PPase_Raf/(pow(KMOLE, 2)*(0.0016611295681063099*B_Raf_active_cyto + PPase_Raf_Km));
            dwdp[19] = -2.7593514420370556e-6*MEK_active_cyto*PP2A_cyto*kcat_PPase_mek/(pow(KMOLE, 2)*(0.0016611295681063099*MEK_active_cyto + PPase_mek_Km));
            dwdp[20] = -2.7593514420370556e-6*MAPK_active_cyto*PTP_cyto*kcat_PTP/(pow(KMOLE, 2)*(0.0016611295681063099*MAPK_active_cyto + PTP_Km));
            dwdp[21] = -2.7593514420370556e-6*MAPK_active_cyto*PTP_PKA_cyto*kcat_PTP_PKA/(pow(KMOLE, 2)*(0.0016611295681063099*MAPK_active_cyto + PTP_PKA_Km));
            dwdp[22] = -2.8973190141389085e-7*B_Raf_active_cyto*MEK_cyto/(pow(KMOLE, 2)*(0.0016611295681063099*MEK_cyto + Raf_activates_MEK_Km));
            dwdp[24] = -(-0.0016611295681063099*GRK_bg_cyto*bg_binds_GRK_Kr_bg_binds_GRK + 2.7593514420370556e-6*GRK_cyto*bg_binds_GRK_Kf_bg_binds_GRK*bg_cyto)/pow(KMOLE, 2);
            dwdp[25] = -2.2074811536296445e-5*PDE_high_km_cyto*cAMP_cyto/(pow(KMOLE, 2)*(0.0016611295681063099*cAMP_cyto + highKM_PDE_Km));
            dwdp[28] = -5.5187028840741121e-5*PDE4_P_cyto*cAMP_cyto/(pow(KMOLE, 2)*(0.0016611295681063099*cAMP_cyto + pde4_p_Km_pde4_p));
            dwdp[29] = -1.379675721018528e-5*PDE4_P_cyto*PP_PDE_cyto/(pow(KMOLE, 2)*(0.0016611295681063099*PDE4_P_cyto + pp2a_4_Km_pp2a_4));
            dwdp[30] = -2.7593514420370556e-6*PTP_PKA_cyto*PTP_PP_cyto*kcat_pp_ptp_pp_ptp/(pow(KMOLE, 2)*(0.0016611295681063099*PTP_PKA_cyto + pp_ptp_Km));
            dwdp[31] = -(2.7593514420370556e-6*G_GDP_cyto*bg_cyto*trimer_Kf_trimer - 0.0016611295681063099*G_protein_cyto*trimer_Kr_trimer)/pow(KMOLE, 2);
            break;
        case 1:
            dwdp[15] = 2.7593514420370556e-6*PKA_cyto*PTP_cyto/(KMOLE*(PKA_P_PTP_Km + 0.0016611295681063099*PTP_cyto));
            break;
        case 2:
            dwdp[16] = 2.7593514420370556e-6*B_Raf_cyto*PKA_cyto/(KMOLE*(0.0016611295681063099*B_Raf_cyto + PKA_activates_Raf_Km));
            break;
        case 3:
            dwdp[17] = 2.7593514420370556e-6*MAPK_active_cyto*PP2A_cyto/(KMOLE*(0.0016611295681063099*MAPK_active_cyto + PPase_MAPK_Km));
            break;
        case 4:
            dwdp[18] = 2.7593514420370556e-6*B_Raf_active_cyto*PP2A_cyto/(KMOLE*(0.0016611295681063099*B_Raf_active_cyto + PPase_Raf_Km));
            break;
        case 5:
            dwdp[19] = 2.7593514420370556e-6*MEK_active_cyto*PP2A_cyto/(KMOLE*(0.0016611295681063099*MEK_active_cyto + PPase_mek_Km));
            break;
        case 6:
            dwdp[20] = 2.7593514420370556e-6*MAPK_active_cyto*PTP_cyto/(KMOLE*(0.0016611295681063099*MAPK_active_cyto + PTP_Km));
            break;
        case 7:
            dwdp[21] = 2.7593514420370556e-6*MAPK_active_cyto*PTP_PKA_cyto/(KMOLE*(0.0016611295681063099*MAPK_active_cyto + PTP_PKA_Km));
            break;
        case 8:
            dwdp[30] = 2.7593514420370556e-6*PTP_PKA_cyto*PTP_PP_cyto/(KMOLE*(0.0016611295681063099*PTP_PKA_cyto + pp_ptp_Km));
            break;
        case 9:
            dwdp[0] = 2.7593514420370556e-6*c2_R2C2_cyto*cAMP_cyto/KMOLE;
            break;
        case 10:
            dwdp[0] = -0.0016611295681063099*c3_R2C2_cyto/KMOLE;
            break;
        case 11:
            dwdp[1] = 2.7593514420370556e-6*c3_R2C2_cyto*cAMP_cyto/KMOLE;
            break;
        case 12:
            dwdp[1] = -0.0016611295681063099*PKA_cyto/KMOLE;
            break;
        case 14:
            dwdp[2] = 0.00033222591362126199*AC_cyto_mem*G_a_s_cyto;
            break;
        case 15:
            dwdp[2] = -0.20000000000000001*AC_active_cyto_mem;
            break;
        case 17:
            dwdp[3] = -0.002823920265780727*AC_active_cyto_mem*ATP_cyto/pow(AC_active_Km_AC_active + 0.0016611295681063099*ATP_cyto, 2);
            break;
        case 19:
            dwdp[4] = -6.6445182724252401e-5*AC_cyto_mem*ATP_cyto/pow(AC_basal_Km_AC_basal + 0.0016611295681063099*ATP_cyto, 2);
            break;
        case 20:
            dwdp[5] = 2.7593514420370556e-6*R2C2_cyto*cAMP_cyto/KMOLE;
            break;
        case 21:
            dwdp[5] = -0.0016611295681063099*c_R2C2_cyto/KMOLE;
            break;
        case 22:
            dwdp[6] = 2.7593514420370556e-6*cAMP_cyto*c_R2C2_cyto/KMOLE;
            break;
        case 23:
            dwdp[6] = -0.0016611295681063099*c2_R2C2_cyto/KMOLE;
            break;
        case 25:
            dwdp[7] = -3.4551495016611247e-5*GRK_cyto*iso_BAR_cyto_mem/pow(GRK_Km_grk + iso_BAR_cyto_mem, 2);
            break;
        case 27:
            dwdp[8] = -0.00044518272425249112*GRK_bg_cyto*iso_BAR_cyto_mem/pow(GRK_bg_Km_GRK_bg + iso_BAR_cyto_mem, 2);
            break;
        case 28:
            dwdp[9] = 0.0016611295681063099*G_a_s_cyto/KMOLE;
            break;
        case 29:
            dwdp[9] = -0.0016611295681063099*G_GDP_cyto/KMOLE;
            break;
        case 31:
            dwdp[10] = 0.00033222591362126199*BAR_cyto_mem*G_protein_cyto;
            break;
        case 32:
            dwdp[10] = -0.20000000000000001*BAR_G_cyto_mem;
            break;
        case 34:
            dwdp[11] = 0.00033222591362126199*G_protein_cyto*iso_BAR_cyto_mem;
            break;
        case 35:
            dwdp[11] = -0.20000000000000001*iso_BAR_G_cyto_mem;
            break;
        case 36:
            dwdp[12] = -4.1390271630555837e-7*MAPK_cyto*MEK_active_cyto/(KMOLE*pow(0.0016611295681063099*MAPK_cyto + MEK_activates_MAPK_Km, 2));
            break;
        case 37:
            dwdp[13] = -2.2074811536296445e-5*PDE4_cyto*cAMP_cyto/(KMOLE*pow(PDE4_Km_PDE4 + 0.0016611295681063099*cAMP_cyto, 2));
            break;
        case 38:
            dwdp[14] = -2.759351442037056e-5*PDE4_cyto*PKA_cyto/(KMOLE*pow(0.0016611295681063099*PDE4_cyto + PKA_P_PDE_Km, 2));
            break;
        case 39:
            dwdp[15] = -2.7593514420370556e-6*PKA_cyto*PTP_cyto*kcat_PKA_P_PTP/(KMOLE*pow(PKA_P_PTP_Km + 0.0016611295681063099*PTP_cyto, 2));
            break;
        case 40:
            dwdp[16] = -2.7593514420370556e-6*B_Raf_cyto*PKA_cyto*kcat_PKA_activates_Raf/(KMOLE*pow(0.0016611295681063099*B_Raf_cyto + PKA_activates_Raf_Km, 2));
            break;
        case 41:
            dwdp[17] = -2.7593514420370556e-6*MAPK_active_cyto*PP2A_cyto*kcat_PPase_MAPK/(KMOLE*pow(0.0016611295681063099*MAPK_active_cyto + PPase_MAPK_Km, 2));
            break;
        case 42:
            dwdp[18] = -2.7593514420370556e-6*B_Raf_active_cyto*PP2A_cyto*kcat_PPase_Raf/(KMOLE*pow(0.0016611295681063099*B_Raf_active_cyto + PPase_Raf_Km, 2));
            break;
        case 43:
            dwdp[19] = -2.7593514420370556e-6*MEK_active_cyto*PP2A_cyto*kcat_PPase_mek/(KMOLE*pow(0.0016611295681063099*MEK_active_cyto + PPase_mek_Km, 2));
            break;
        case 44:
            dwdp[20] = -2.7593514420370556e-6*MAPK_active_cyto*PTP_cyto*kcat_PTP/(KMOLE*pow(0.0016611295681063099*MAPK_active_cyto + PTP_Km, 2));
            break;
        case 45:
            dwdp[21] = -2.7593514420370556e-6*MAPK_active_cyto*PTP_PKA_cyto*kcat_PTP_PKA/(KMOLE*pow(0.0016611295681063099*MAPK_active_cyto + PTP_PKA_Km, 2));
            break;
        case 46:
            dwdp[22] = -2.8973190141389085e-7*B_Raf_active_cyto*MEK_cyto/(KMOLE*pow(0.0016611295681063099*MEK_cyto + Raf_activates_MEK_Km, 2));
            break;
        case 48:
            dwdp[23] = 0.20000000000000001*iso_BAR_G_cyto_mem;
            break;
        case 49:
            dwdp[23] = -5.5187028840741116e-7*G_a_s_cyto*bg_cyto*iso_BAR_cyto_mem;
            break;
        case 50:
            dwdp[24] = 2.7593514420370556e-6*GRK_cyto*bg_cyto/KMOLE;
            break;
        case 51:
            dwdp[24] = -0.0016611295681063099*GRK_bg_cyto/KMOLE;
            break;
        case 52:
            dwdp[25] = -2.2074811536296445e-5*PDE_high_km_cyto*cAMP_cyto/(KMOLE*pow(0.0016611295681063099*cAMP_cyto + highKM_PDE_Km, 2));
            break;
        case 54:
            dwdp[26] = 0.00033222591362126199*BAR_cyto_mem*iso_extra;
            break;
        case 55:
            dwdp[26] = -0.20000000000000001*iso_BAR_cyto_mem;
            break;
        case 57:
            dwdp[27] = 0.00033222591362126199*BAR_G_cyto_mem*iso_extra;
            break;
        case 58:
            dwdp[27] = -0.20000000000000001*iso_BAR_G_cyto_mem;
            break;
        case 59:
            dwdp[28] = -5.5187028840741121e-5*PDE4_P_cyto*cAMP_cyto/(KMOLE*pow(0.0016611295681063099*cAMP_cyto + pde4_p_Km_pde4_p, 2));
            break;
        case 60:
            dwdp[29] = -1.379675721018528e-5*PDE4_P_cyto*PP_PDE_cyto/(KMOLE*pow(0.0016611295681063099*PDE4_P_cyto + pp2a_4_Km_pp2a_4, 2));
            break;
        case 61:
            dwdp[30] = -2.7593514420370556e-6*PTP_PKA_cyto*PTP_PP_cyto*kcat_pp_ptp_pp_ptp/(KMOLE*pow(0.0016611295681063099*PTP_PKA_cyto + pp_ptp_Km, 2));
            break;
        case 62:
            dwdp[31] = 2.7593514420370556e-6*G_GDP_cyto*bg_cyto/KMOLE;
            break;
        case 63:
            dwdp[31] = -0.0016611295681063099*G_protein_cyto/KMOLE;
            break;
    }
}