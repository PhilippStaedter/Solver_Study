#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_jiang1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*GAP;
            dwdp[1] = -1.0*GLC + 0.001*GLCflow_Glc_F;
            dwdp[2] = 1.0*LAC;
            break;
        case 1:
            dwdp[1] = 0.001*flow;
            break;
        case 2:
            dwdp[3] = -1.0*pow(ADP_cyt, 2);
            break;
        case 3:
            dwdp[3] = 1.0*AMP*ATP_cyt;
            break;
        case 4:
            dwdp[4] = -1.0*ATP_cyt*GLC*v1_V1/(pow(ATP_cyt + v1_K1ATP, 2)*(GLC + v1_K1GLC));
            break;
        case 5:
            dwdp[4] = -1.0*ATP_cyt*GLC*v1_V1/((ATP_cyt + v1_K1ATP)*pow(GLC + v1_K1GLC, 2));
            break;
        case 6:
            dwdp[4] = 1.0*ATP_cyt*GLC/((ATP_cyt + v1_K1ATP)*(GLC + v1_K1GLC));
            break;
        case 7:
            dwdp[5] = 1.0*Acetyl_CoA*OXA*v10_V/(Acetyl_CoA*OXA + Acetyl_CoA*v10_Kb + OXA*v10_Ka + v10_Kia*v10_Kib);
            break;
        case 8:
            dwdp[5] = -1.0*Acetyl_CoA*OXA*v10_Kia*v10_V*v10_v10_CS/pow(Acetyl_CoA*OXA + Acetyl_CoA*v10_Kb + OXA*v10_Ka + v10_Kia*v10_Kib, 2);
            break;
        case 9:
            dwdp[5] = -1.0*Acetyl_CoA*OXA*v10_Kib*v10_V*v10_v10_CS/pow(Acetyl_CoA*OXA + Acetyl_CoA*v10_Kb + OXA*v10_Ka + v10_Kia*v10_Kib, 2);
            break;
        case 10:
            dwdp[5] = -1.0*pow(Acetyl_CoA, 2)*OXA*v10_V*v10_v10_CS/pow(Acetyl_CoA*OXA + Acetyl_CoA*v10_Kb + OXA*v10_Ka + v10_Kia*v10_Kib, 2);
            break;
        case 11:
            dwdp[5] = -1.0*Acetyl_CoA*pow(OXA, 2)*v10_V*v10_v10_CS/pow(Acetyl_CoA*OXA + Acetyl_CoA*v10_Kb + OXA*v10_Ka + v10_Kia*v10_Kib, 2);
            break;
        case 12:
            dwdp[5] = 1.0*Acetyl_CoA*OXA*v10_v10_CS/(Acetyl_CoA*OXA + Acetyl_CoA*v10_Kb + OXA*v10_Ka + v10_Kia*v10_Kib);
            break;
        case 13:
            dwdp[6] = 1.0*(Cit*v11_KcF*v11_Kp - IsoCit*v11_KcR*v11_Ks)/(Cit*v11_Kp + IsoCit*v11_Ks + v11_Kp*v11_Ks);
            break;
        case 14:
            dwdp[6] = -1.0*IsoCit*v11_Ks*v11_v11_ACO/(Cit*v11_Kp + IsoCit*v11_Ks + v11_Kp*v11_Ks);
            break;
        case 15:
            dwdp[6] = 1.0*Cit*v11_Kp*v11_v11_ACO/(Cit*v11_Kp + IsoCit*v11_Ks + v11_Kp*v11_Ks);
            break;
        case 16:
            dwdp[6] = 1.0*Cit*v11_KcF*v11_v11_ACO/(Cit*v11_Kp + IsoCit*v11_Ks + v11_Kp*v11_Ks) + 1.0*v11_v11_ACO*(-Cit - v11_Ks)*(Cit*v11_KcF*v11_Kp - IsoCit*v11_KcR*v11_Ks)/pow(Cit*v11_Kp + IsoCit*v11_Ks + v11_Kp*v11_Ks, 2);
            break;
        case 17:
            dwdp[6] = -1.0*IsoCit*v11_KcR*v11_v11_ACO/(Cit*v11_Kp + IsoCit*v11_Ks + v11_Kp*v11_Ks) + 1.0*v11_v11_ACO*(-IsoCit - v11_Kp)*(Cit*v11_KcF*v11_Kp - IsoCit*v11_KcR*v11_Ks)/pow(Cit*v11_Kp + IsoCit*v11_Ks + v11_Kp*v11_Ks, 2);
            break;
        case 18:
            dwdp[7] = 1.0*v12_KcF*(ADP*IsoCit*v12_b + pow(IsoCit, 2))/(ADP*IsoCit*v12_e + ADP*v12_d + pow(IsoCit, 2) + IsoCit*v12_c + v12_f);
            break;
        case 19:
            dwdp[7] = -1.0*v12_KcF*v12_v12_IDHa*(ADP*IsoCit*v12_b + pow(IsoCit, 2))/pow(ADP*IsoCit*v12_e + ADP*v12_d + pow(IsoCit, 2) + IsoCit*v12_c + v12_f, 2);
            break;
        case 20:
            dwdp[7] = -1.0*ADP*IsoCit*v12_KcF*v12_v12_IDHa*(ADP*IsoCit*v12_b + pow(IsoCit, 2))/pow(ADP*IsoCit*v12_e + ADP*v12_d + pow(IsoCit, 2) + IsoCit*v12_c + v12_f, 2);
            break;
        case 21:
            dwdp[7] = -1.0*ADP*v12_KcF*v12_v12_IDHa*(ADP*IsoCit*v12_b + pow(IsoCit, 2))/pow(ADP*IsoCit*v12_e + ADP*v12_d + pow(IsoCit, 2) + IsoCit*v12_c + v12_f, 2);
            break;
        case 22:
            dwdp[7] = -1.0*IsoCit*v12_KcF*v12_v12_IDHa*(ADP*IsoCit*v12_b + pow(IsoCit, 2))/pow(ADP*IsoCit*v12_e + ADP*v12_d + pow(IsoCit, 2) + IsoCit*v12_c + v12_f, 2);
            break;
        case 23:
            dwdp[7] = 1.0*ADP*IsoCit*v12_KcF*v12_v12_IDHa/(ADP*IsoCit*v12_e + ADP*v12_d + pow(IsoCit, 2) + IsoCit*v12_c + v12_f);
            break;
        case 24:
            dwdp[7] = 1.0*v12_v12_IDHa*(ADP*IsoCit*v12_b + pow(IsoCit, 2))/(ADP*IsoCit*v12_e + ADP*v12_d + pow(IsoCit, 2) + IsoCit*v12_c + v12_f);
            break;
        case 25:
            dwdp[8] = 1.0*CoA*NAD_p*OG*v14_KcF/(CoA*NADH*OG*v14_KmC/v14_Kir + CoA*NAD_p*OG + CoA*NAD_p*v14_KmA + CoA*OG*v14_KmC + NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR) + NAD_p*OG*SCoA*v14_KmB/v14_Kiq + NAD_p*OG*v14_KmB);
            break;
        case 26:
            dwdp[8] = 1.0*CoA*NAD_p*OG*v14_v14_OGDC/(CoA*NADH*OG*v14_KmC/v14_Kir + CoA*NAD_p*OG + CoA*NAD_p*v14_KmA + CoA*OG*v14_KmC + NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR) + NAD_p*OG*SCoA*v14_KmB/v14_Kiq + NAD_p*OG*v14_KmB);
            break;
        case 27:
            dwdp[8] = 1.0*pow(CoA, 2)*NADH*NAD_p*pow(OG, 2)*v14_KcF*v14_KmC*v14_v14_OGDC/(pow(v14_Kir, 2)*pow(CoA*NADH*OG*v14_KmC/v14_Kir + CoA*NAD_p*OG + CoA*NAD_p*v14_KmA + CoA*OG*v14_KmC + NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR) + NAD_p*OG*SCoA*v14_KmB/v14_Kiq + NAD_p*OG*v14_KmB, 2));
            break;
        case 28:
            dwdp[8] = 1.0*CoA*NAD_p*OG*v14_KcF*v14_v14_OGDC*(NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*pow(v14_Kiq, 2)*v14_KmR) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*pow(v14_Kiq, 2)*v14_KmR) + NAD_p*OG*SCoA*v14_KmB/pow(v14_Kiq, 2))/pow(CoA*NADH*OG*v14_KmC/v14_Kir + CoA*NAD_p*OG + CoA*NAD_p*v14_KmA + CoA*OG*v14_KmC + NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR) + NAD_p*OG*SCoA*v14_KmB/v14_Kiq + NAD_p*OG*v14_KmB, 2);
            break;
        case 29:
            dwdp[8] = 1.0*CoA*NAD_p*OG*v14_KcF*v14_v14_OGDC*(NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*pow(v14_Kip, 2)*v14_Kiq*v14_KmR) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(pow(v14_Kip, 2)*v14_Kiq*v14_KmR))/pow(CoA*NADH*OG*v14_KmC/v14_Kir + CoA*NAD_p*OG + CoA*NAD_p*v14_KmA + CoA*OG*v14_KmC + NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR) + NAD_p*OG*SCoA*v14_KmB/v14_Kiq + NAD_p*OG*v14_KmB, 2);
            break;
        case 30:
            dwdp[8] = 1.0*CoA*NAD_p*OG*v14_KcF*v14_v14_OGDC*(-NADH*OG*SCoA*v14_Kib*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) - NADH*SCoA*v14_Kib*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR))/pow(CoA*NADH*OG*v14_KmC/v14_Kir + CoA*NAD_p*OG + CoA*NAD_p*v14_KmA + CoA*OG*v14_KmC + NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR) + NAD_p*OG*SCoA*v14_KmB/v14_Kiq + NAD_p*OG*v14_KmB, 2);
            break;
        case 31:
            dwdp[8] = 1.0*CoA*NAD_p*OG*v14_KcF*v14_v14_OGDC*(-NADH*OG*SCoA*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) - NADH*SCoA*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR))/pow(CoA*NADH*OG*v14_KmC/v14_Kir + CoA*NAD_p*OG + CoA*NAD_p*v14_KmA + CoA*OG*v14_KmC + NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR) + NAD_p*OG*SCoA*v14_KmB/v14_Kiq + NAD_p*OG*v14_KmB, 2);
            break;
        case 32:
            dwdp[8] = 1.0*CoA*NADH*NAD_p*pow(OG, 2)*SCoA*v14_KcF*v14_Kib*v14_Kic*v14_KmA*v14_KmP*v14_v14_OGDC/(pow(v14_Kia, 2)*v14_Kip*v14_Kiq*v14_KmR*pow(CoA*NADH*OG*v14_KmC/v14_Kir + CoA*NAD_p*OG + CoA*NAD_p*v14_KmA + CoA*OG*v14_KmC + NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR) + NAD_p*OG*SCoA*v14_KmB/v14_Kiq + NAD_p*OG*v14_KmB, 2));
            break;
        case 33:
            dwdp[8] = 1.0*CoA*NAD_p*OG*v14_KcF*v14_v14_OGDC*(NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*pow(v14_KmR, 2)) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*pow(v14_KmR, 2)))/pow(CoA*NADH*OG*v14_KmC/v14_Kir + CoA*NAD_p*OG + CoA*NAD_p*v14_KmA + CoA*OG*v14_KmC + NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR) + NAD_p*OG*SCoA*v14_KmB/v14_Kiq + NAD_p*OG*v14_KmB, 2);
            break;
        case 34:
            dwdp[8] = 1.0*CoA*NAD_p*OG*v14_KcF*v14_v14_OGDC*(-NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) - NADH*SCoA*v14_Kib*v14_Kic*v14_KmA/(v14_Kip*v14_Kiq*v14_KmR))/pow(CoA*NADH*OG*v14_KmC/v14_Kir + CoA*NAD_p*OG + CoA*NAD_p*v14_KmA + CoA*OG*v14_KmC + NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR) + NAD_p*OG*SCoA*v14_KmB/v14_Kiq + NAD_p*OG*v14_KmB, 2);
            break;
        case 35:
            dwdp[8] = 1.0*CoA*NAD_p*OG*v14_KcF*v14_v14_OGDC*(-CoA*NADH*OG/v14_Kir - CoA*OG)/pow(CoA*NADH*OG*v14_KmC/v14_Kir + CoA*NAD_p*OG + CoA*NAD_p*v14_KmA + CoA*OG*v14_KmC + NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR) + NAD_p*OG*SCoA*v14_KmB/v14_Kiq + NAD_p*OG*v14_KmB, 2);
            break;
        case 36:
            dwdp[8] = 1.0*CoA*NAD_p*OG*v14_KcF*v14_v14_OGDC*(-NAD_p*OG*SCoA/v14_Kiq - NAD_p*OG)/pow(CoA*NADH*OG*v14_KmC/v14_Kir + CoA*NAD_p*OG + CoA*NAD_p*v14_KmA + CoA*OG*v14_KmC + NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR) + NAD_p*OG*SCoA*v14_KmB/v14_Kiq + NAD_p*OG*v14_KmB, 2);
            break;
        case 37:
            dwdp[8] = 1.0*CoA*NAD_p*OG*v14_KcF*v14_v14_OGDC*(-CoA*NAD_p - NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) - NADH*SCoA*v14_Kib*v14_Kic*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR))/pow(CoA*NADH*OG*v14_KmC/v14_Kir + CoA*NAD_p*OG + CoA*NAD_p*v14_KmA + CoA*OG*v14_KmC + NADH*OG*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kia*v14_Kip*v14_Kiq*v14_KmR) + NADH*SCoA*v14_Kib*v14_Kic*v14_KmA*v14_KmP/(v14_Kip*v14_Kiq*v14_KmR) + NAD_p*OG*SCoA*v14_KmB/v14_Kiq + NAD_p*OG*v14_KmB, 2);
            break;
        case 38:
            dwdp[9] = 1.0*(v15_Kc1 + v15_Kc2*(Suc*v15_Kip*v15_KmC/v15_KmC2 + amici::pi/v15_KmC2))*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)/(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB);
            break;
        case 39:
            dwdp[9] = 1.0*v15_v15_SCS*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)*(Suc*v15_Kip*v15_KmC/v15_KmC2 + amici::pi/v15_KmC2)/(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB);
            break;
        case 40:
            dwdp[9] = 1.0*v15_v15_SCS*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)/(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB);
            break;
        case 41:
            dwdp[9] = 1.0*(v15_Kc1*v15_v15_SCS + v15_Kc2*v15_v15_SCS*(Suc*v15_Kip*v15_KmC/v15_KmC2 + amici::pi/v15_KmC2))*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)*(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*pow(v15_Kir, 2)) + amici::pi*CoA*GDP*v15_KmB/pow(v15_Kir, 2) + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*pow(v15_Kir, 2)*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(pow(v15_Kir, 2)*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*pow(v15_Kir, 2)*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(pow(v15_Kir, 2)*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(pow(v15_Kir, 2)*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*pow(v15_Kir, 2)) + amici::pi*CoA*v15_Kia*v15_KmB/pow(v15_Kir, 2))/pow(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB, 2);
            break;
        case 42:
            dwdp[9] = 1.0*(v15_Kc1*v15_v15_SCS + v15_Kc2*v15_v15_SCS*(Suc*v15_Kip*v15_KmC/v15_KmC2 + amici::pi/v15_KmC2))*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)*(GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*pow(v15_Kiq, 2)) + amici::pi*GTP*SCoA*v15_KmA/pow(v15_Kiq, 2) + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*pow(v15_Kiq, 2)) + amici::pi*GTP*v15_Kia*v15_KmB/pow(v15_Kiq, 2))/pow(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB, 2);
            break;
        case 43:
            dwdp[9] = 1.0*Suc*v15_Kc2*v15_KmC*v15_v15_SCS*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)/(v15_KmC2*(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB)) + 1.0*(v15_Kc1*v15_v15_SCS + v15_Kc2*v15_v15_SCS*(Suc*v15_Kip*v15_KmC/v15_KmC2 + amici::pi/v15_KmC2))*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)*(CoA*GDP*Suc*v15_KmB*v15_KmC/(pow(v15_Kip, 2)*v15_Kir) + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(pow(v15_Kip, 2)*v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(pow(v15_Kip, 2)*v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(pow(v15_Kip, 2)*v15_Kir) + GDP*SCoA*Suc*v15_KmC/pow(v15_Kip, 2) + amici::pi*GDP*SCoA*Suc*v15_KmC/(pow(v15_Kip, 2)*v15_KmC2) + GDP*Suc*v15_KmB*v15_KmC/pow(v15_Kip, 2) + GTP*SCoA*Suc*v15_KmA*v15_KmC/(pow(v15_Kip, 2)*v15_Kiq) + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(pow(v15_Kip, 2)*v15_Kiq) + SCoA*Suc*v15_KmA*v15_KmC/pow(v15_Kip, 2) + Suc*v15_Kia*v15_KmB*v15_KmC/pow(v15_Kip, 2))/pow(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB, 2);
            break;
        case 44:
            dwdp[9] = -1.0*CoA*GTP*v15_Kia*v15_KmB*(v15_Kc1*v15_v15_SCS + v15_Kc2*v15_v15_SCS*(Suc*v15_Kip*v15_KmC/v15_KmC2 + amici::pi/v15_KmC2))*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)/(v15_Kir*v15_KmQ*pow(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB, 2));
            break;
        case 46:
            dwdp[9] = 1.0*(v15_Kc1*v15_v15_SCS + v15_Kc2*v15_v15_SCS*(Suc*v15_Kip*v15_KmC/v15_KmC2 + amici::pi/v15_KmC2))*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)*(-CoA*GTP*pow(Suc, 2)*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) - amici::pi*CoA*GTP*Suc*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) - CoA*GTP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) - CoA*GTP*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) - amici::pi*CoA*GTP*v15_KmB/(v15_Kir*v15_KmQ) - CoA*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) - amici::pi*CoA*v15_KmB/v15_Kir - GTP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) - amici::pi*GTP*v15_KmB/v15_Kiq - Suc*v15_KmB*v15_KmC/v15_Kip - amici::pi*v15_KmB)/pow(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB, 2);
            break;
        case 47:
            dwdp[9] = 1.0*CoA*GTP*Suc*(v15_Kc1*v15_v15_SCS + v15_Kc2*v15_v15_SCS*(Suc*v15_Kip*v15_KmC/v15_KmC2 + amici::pi/v15_KmC2))/(pow(v15_Keq, 2)*(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB));
            break;
        case 48:
            dwdp[9] = 1.0*(v15_Kc1*v15_v15_SCS + v15_Kc2*v15_v15_SCS*(Suc*v15_Kip*v15_KmC/v15_KmC2 + amici::pi/v15_KmC2))*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)*(CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*pow(v15_KmP2, 2)*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*pow(v15_KmP2, 2)*v15_KmQ))/pow(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB, 2);
            break;
        case 49:
            dwdp[9] = 1.0*v15_Kc2*v15_v15_SCS*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)*(-Suc*v15_Kip*v15_KmC/pow(v15_KmC2, 2) - amici::pi/pow(v15_KmC2, 2))/(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB) + 1.0*(v15_Kc1*v15_v15_SCS + v15_Kc2*v15_v15_SCS*(Suc*v15_Kip*v15_KmC/v15_KmC2 + amici::pi/v15_KmC2))*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)*(amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*pow(v15_KmC2, 2)) + pow(amici::pi, 2)*GDP*SCoA/pow(v15_KmC2, 2))/pow(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB, 2);
            break;
        case 50:
            dwdp[9] = 1.0*(v15_Kc1*v15_v15_SCS + v15_Kc2*v15_v15_SCS*(Suc*v15_Kip*v15_KmC/v15_KmC2 + amici::pi/v15_KmC2))*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)*(CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*pow(v15_KmQ, 2)) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*pow(v15_KmQ, 2)) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*pow(v15_KmQ, 2)) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*pow(v15_KmQ, 2)) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*pow(v15_KmQ, 2)))/pow(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB, 2);
            break;
        case 52:
            dwdp[9] = 1.0*Suc*v15_Kc2*v15_Kip*v15_v15_SCS*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)/(v15_KmC2*(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB)) + 1.0*(v15_Kc1*v15_v15_SCS + v15_Kc2*v15_v15_SCS*(Suc*v15_Kip*v15_KmC/v15_KmC2 + amici::pi/v15_KmC2))*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)*(-CoA*GDP*Suc*v15_KmB/(v15_Kip*v15_Kir) - CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) - CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kip*v15_Kir*v15_KmQ) - CoA*Suc*v15_Kia*v15_KmB/(v15_Kip*v15_Kir) - GDP*SCoA*Suc/v15_Kip - amici::pi*GDP*SCoA*Suc/(v15_Kip*v15_KmC2) - GDP*SCoA - GDP*Suc*v15_KmB/v15_Kip - GTP*SCoA*Suc*v15_KmA/(v15_Kip*v15_Kiq) - GTP*Suc*v15_Kia*v15_KmB/(v15_Kip*v15_Kiq) - SCoA*Suc*v15_KmA/v15_Kip - Suc*v15_Kia*v15_KmB/v15_Kip)/pow(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB, 2);
            break;
        case 53:
            dwdp[9] = 1.0*(v15_Kc1*v15_v15_SCS + v15_Kc2*v15_v15_SCS*(Suc*v15_Kip*v15_KmC/v15_KmC2 + amici::pi/v15_KmC2))*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)*(-CoA*GDP*Suc*v15_KmC/(v15_Kip*v15_Kir) - amici::pi*CoA*GDP/v15_Kir - CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) - amici::pi*CoA*GTP*Suc*v15_Kia/(v15_Kir*v15_KmP2*v15_KmQ) - CoA*GTP*Suc*v15_Kia*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) - CoA*GTP*v15_Kia*v15_Kic/(v15_Kir*v15_KmQ) - amici::pi*CoA*GTP*v15_Kia/(v15_Kir*v15_KmQ) - CoA*Suc*v15_Kia*v15_KmC/(v15_Kip*v15_Kir) - amici::pi*CoA*v15_Kia/v15_Kir - GDP*Suc*v15_KmC/v15_Kip - amici::pi*GDP - GTP*Suc*v15_Kia*v15_KmC/(v15_Kip*v15_Kiq) - amici::pi*GTP*v15_Kia/v15_Kiq - Suc*v15_Kia*v15_KmC/v15_Kip - amici::pi*v15_Kia)/pow(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB, 2);
            break;
        case 54:
            dwdp[9] = 1.0*(v15_Kc1*v15_v15_SCS + v15_Kc2*v15_v15_SCS*(Suc*v15_Kip*v15_KmC/v15_KmC2 + amici::pi/v15_KmC2))*(-CoA*GTP*Suc/v15_Keq + amici::pi*GDP*SCoA)*(-GTP*SCoA*Suc*v15_KmC/(v15_Kip*v15_Kiq) - amici::pi*GTP*SCoA/v15_Kiq - SCoA*Suc*v15_KmC/v15_Kip - amici::pi*SCoA)/pow(CoA*GDP*Suc*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*GDP*v15_KmB/v15_Kir + CoA*GTP*pow(Suc, 2)*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmP2*v15_KmQ) + amici::pi*CoA*GTP*Suc*v15_Kia*v15_KmB/(v15_Kir*v15_KmP2*v15_KmQ) + CoA*GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir*v15_KmQ) + CoA*GTP*v15_Kia*v15_Kic*v15_KmB/(v15_Kir*v15_KmQ) + amici::pi*CoA*GTP*v15_Kia*v15_KmB/(v15_Kir*v15_KmQ) + CoA*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kir) + amici::pi*CoA*v15_Kia*v15_KmB/v15_Kir + GDP*SCoA*Suc*v15_KmC/v15_Kip + amici::pi*GDP*SCoA*Suc*v15_KmC/(v15_Kip*v15_KmC2) + GDP*SCoA*v15_KmC + amici::pi*GDP*SCoA + pow(amici::pi, 2)*GDP*SCoA/v15_KmC2 + GDP*Suc*v15_KmB*v15_KmC/v15_Kip + amici::pi*GDP*v15_KmB + GTP*SCoA*Suc*v15_KmA*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*SCoA*v15_KmA/v15_Kiq + GTP*Suc*v15_Kia*v15_KmB*v15_KmC/(v15_Kip*v15_Kiq) + amici::pi*GTP*v15_Kia*v15_KmB/v15_Kiq + SCoA*Suc*v15_KmA*v15_KmC/v15_Kip + amici::pi*SCoA*v15_KmA + Suc*v15_Kia*v15_KmB*v15_KmC/v15_Kip + amici::pi*v15_Kia*v15_KmB, 2);
            break;
        case 55:
            dwdp[10] = 1.0*v16_KcF*v16_KcR*(-Fum*QH2/v16_Keq + Q*Suc)/(Fum*QH2*v16_KcF/v16_Keq + Fum*Suc*v16_KcF*v16_KmP2/(v16_Keq*v16_KiS1) + Fum*v16_KcF*v16_KmP2/v16_Keq + Q*QH2*v16_KcR*v16_KmS1/v16_KiP2 + Q*Suc*v16_KcR + Q*v16_KcR*v16_KmS1 + QH2*v16_KcF*v16_KmP1/v16_Keq + Suc*v16_KcR*v16_KmS2);
            break;
        case 56:
            dwdp[10] = 1.0*v16_KcF*v16_KcR*v16_v16_SDH*(-Fum*QH2/v16_Keq + Q*Suc)*(-Q*QH2*v16_KmS1/v16_KiP2 - Q*Suc - Q*v16_KmS1 - Suc*v16_KmS2)/pow(Fum*QH2*v16_KcF/v16_Keq + Fum*Suc*v16_KcF*v16_KmP2/(v16_Keq*v16_KiS1) + Fum*v16_KcF*v16_KmP2/v16_Keq + Q*QH2*v16_KcR*v16_KmS1/v16_KiP2 + Q*Suc*v16_KcR + Q*v16_KcR*v16_KmS1 + QH2*v16_KcF*v16_KmP1/v16_Keq + Suc*v16_KcR*v16_KmS2, 2) + 1.0*v16_KcF*v16_v16_SDH*(-Fum*QH2/v16_Keq + Q*Suc)/(Fum*QH2*v16_KcF/v16_Keq + Fum*Suc*v16_KcF*v16_KmP2/(v16_Keq*v16_KiS1) + Fum*v16_KcF*v16_KmP2/v16_Keq + Q*QH2*v16_KcR*v16_KmS1/v16_KiP2 + Q*Suc*v16_KcR + Q*v16_KcR*v16_KmS1 + QH2*v16_KcF*v16_KmP1/v16_Keq + Suc*v16_KcR*v16_KmS2);
            break;
        case 57:
            dwdp[10] = 1.0*v16_KcF*v16_KcR*v16_v16_SDH*(-Fum*QH2/v16_Keq + Q*Suc)*(-Fum*QH2/v16_Keq - Fum*Suc*v16_KmP2/(v16_Keq*v16_KiS1) - Fum*v16_KmP2/v16_Keq - QH2*v16_KmP1/v16_Keq)/pow(Fum*QH2*v16_KcF/v16_Keq + Fum*Suc*v16_KcF*v16_KmP2/(v16_Keq*v16_KiS1) + Fum*v16_KcF*v16_KmP2/v16_Keq + Q*QH2*v16_KcR*v16_KmS1/v16_KiP2 + Q*Suc*v16_KcR + Q*v16_KcR*v16_KmS1 + QH2*v16_KcF*v16_KmP1/v16_Keq + Suc*v16_KcR*v16_KmS2, 2) + 1.0*v16_KcR*v16_v16_SDH*(-Fum*QH2/v16_Keq + Q*Suc)/(Fum*QH2*v16_KcF/v16_Keq + Fum*Suc*v16_KcF*v16_KmP2/(v16_Keq*v16_KiS1) + Fum*v16_KcF*v16_KmP2/v16_Keq + Q*QH2*v16_KcR*v16_KmS1/v16_KiP2 + Q*Suc*v16_KcR + Q*v16_KcR*v16_KmS1 + QH2*v16_KcF*v16_KmP1/v16_Keq + Suc*v16_KcR*v16_KmS2);
            break;
        case 58:
            dwdp[10] = 1.0*Fum*QH2*v16_KcF*v16_KcR*v16_v16_SDH/(pow(v16_Keq, 2)*(Fum*QH2*v16_KcF/v16_Keq + Fum*Suc*v16_KcF*v16_KmP2/(v16_Keq*v16_KiS1) + Fum*v16_KcF*v16_KmP2/v16_Keq + Q*QH2*v16_KcR*v16_KmS1/v16_KiP2 + Q*Suc*v16_KcR + Q*v16_KcR*v16_KmS1 + QH2*v16_KcF*v16_KmP1/v16_Keq + Suc*v16_KcR*v16_KmS2)) + 1.0*v16_KcF*v16_KcR*v16_v16_SDH*(-Fum*QH2/v16_Keq + Q*Suc)*(Fum*QH2*v16_KcF/pow(v16_Keq, 2) + Fum*Suc*v16_KcF*v16_KmP2/(pow(v16_Keq, 2)*v16_KiS1) + Fum*v16_KcF*v16_KmP2/pow(v16_Keq, 2) + QH2*v16_KcF*v16_KmP1/pow(v16_Keq, 2))/pow(Fum*QH2*v16_KcF/v16_Keq + Fum*Suc*v16_KcF*v16_KmP2/(v16_Keq*v16_KiS1) + Fum*v16_KcF*v16_KmP2/v16_Keq + Q*QH2*v16_KcR*v16_KmS1/v16_KiP2 + Q*Suc*v16_KcR + Q*v16_KcR*v16_KmS1 + QH2*v16_KcF*v16_KmP1/v16_Keq + Suc*v16_KcR*v16_KmS2, 2);
            break;
        case 59:
            dwdp[10] = 1.0*Q*QH2*v16_KcF*pow(v16_KcR, 2)*v16_KmS1*v16_v16_SDH*(-Fum*QH2/v16_Keq + Q*Suc)/(pow(v16_KiP2, 2)*pow(Fum*QH2*v16_KcF/v16_Keq + Fum*Suc*v16_KcF*v16_KmP2/(v16_Keq*v16_KiS1) + Fum*v16_KcF*v16_KmP2/v16_Keq + Q*QH2*v16_KcR*v16_KmS1/v16_KiP2 + Q*Suc*v16_KcR + Q*v16_KcR*v16_KmS1 + QH2*v16_KcF*v16_KmP1/v16_Keq + Suc*v16_KcR*v16_KmS2, 2));
            break;
        case 60:
            dwdp[10] = 1.0*Fum*Suc*pow(v16_KcF, 2)*v16_KcR*v16_KmP2*v16_v16_SDH*(-Fum*QH2/v16_Keq + Q*Suc)/(v16_Keq*pow(v16_KiS1, 2)*pow(Fum*QH2*v16_KcF/v16_Keq + Fum*Suc*v16_KcF*v16_KmP2/(v16_Keq*v16_KiS1) + Fum*v16_KcF*v16_KmP2/v16_Keq + Q*QH2*v16_KcR*v16_KmS1/v16_KiP2 + Q*Suc*v16_KcR + Q*v16_KcR*v16_KmS1 + QH2*v16_KcF*v16_KmP1/v16_Keq + Suc*v16_KcR*v16_KmS2, 2));
            break;
        case 61:
            dwdp[10] = 1.0*v16_KcF*v16_KcR*v16_v16_SDH*(-Fum*QH2/v16_Keq + Q*Suc)*(-Fum*Suc*v16_KcF/(v16_Keq*v16_KiS1) - Fum*v16_KcF/v16_Keq)/pow(Fum*QH2*v16_KcF/v16_Keq + Fum*Suc*v16_KcF*v16_KmP2/(v16_Keq*v16_KiS1) + Fum*v16_KcF*v16_KmP2/v16_Keq + Q*QH2*v16_KcR*v16_KmS1/v16_KiP2 + Q*Suc*v16_KcR + Q*v16_KcR*v16_KmS1 + QH2*v16_KcF*v16_KmP1/v16_Keq + Suc*v16_KcR*v16_KmS2, 2);
            break;
        case 62:
            dwdp[10] = -1.0*QH2*pow(v16_KcF, 2)*v16_KcR*v16_v16_SDH*(-Fum*QH2/v16_Keq + Q*Suc)/(v16_Keq*pow(Fum*QH2*v16_KcF/v16_Keq + Fum*Suc*v16_KcF*v16_KmP2/(v16_Keq*v16_KiS1) + Fum*v16_KcF*v16_KmP2/v16_Keq + Q*QH2*v16_KcR*v16_KmS1/v16_KiP2 + Q*Suc*v16_KcR + Q*v16_KcR*v16_KmS1 + QH2*v16_KcF*v16_KmP1/v16_Keq + Suc*v16_KcR*v16_KmS2, 2));
            break;
        case 63:
            dwdp[10] = -1.0*Suc*v16_KcF*pow(v16_KcR, 2)*v16_v16_SDH*(-Fum*QH2/v16_Keq + Q*Suc)/pow(Fum*QH2*v16_KcF/v16_Keq + Fum*Suc*v16_KcF*v16_KmP2/(v16_Keq*v16_KiS1) + Fum*v16_KcF*v16_KmP2/v16_Keq + Q*QH2*v16_KcR*v16_KmS1/v16_KiP2 + Q*Suc*v16_KcR + Q*v16_KcR*v16_KmS1 + QH2*v16_KcF*v16_KmP1/v16_Keq + Suc*v16_KcR*v16_KmS2, 2);
            break;
        case 64:
            dwdp[10] = 1.0*v16_KcF*v16_KcR*v16_v16_SDH*(-Fum*QH2/v16_Keq + Q*Suc)*(-Q*QH2*v16_KcR/v16_KiP2 - Q*v16_KcR)/pow(Fum*QH2*v16_KcF/v16_Keq + Fum*Suc*v16_KcF*v16_KmP2/(v16_Keq*v16_KiS1) + Fum*v16_KcF*v16_KmP2/v16_Keq + Q*QH2*v16_KcR*v16_KmS1/v16_KiP2 + Q*Suc*v16_KcR + Q*v16_KcR*v16_KmS1 + QH2*v16_KcF*v16_KmP1/v16_Keq + Suc*v16_KcR*v16_KmS2, 2);
            break;
        case 65:
            dwdp[11] = 1.0*(Fum*v17_KcF*v17_Kp - Mal*v17_KcR*v17_Ks)/(Fum*v17_Kp + Mal*v17_Ks + v17_Kp*v17_Ks);
            break;
        case 66:
            dwdp[11] = -1.0*Mal*v17_Ks*v17_v17_FM/(Fum*v17_Kp + Mal*v17_Ks + v17_Kp*v17_Ks);
            break;
        case 67:
            dwdp[11] = 1.0*Fum*v17_Kp*v17_v17_FM/(Fum*v17_Kp + Mal*v17_Ks + v17_Kp*v17_Ks);
            break;
        case 68:
            dwdp[11] = -1.0*Mal*v17_KcR*v17_v17_FM/(Fum*v17_Kp + Mal*v17_Ks + v17_Kp*v17_Ks) + 1.0*v17_v17_FM*(-Mal - v17_Kp)*(Fum*v17_KcF*v17_Kp - Mal*v17_KcR*v17_Ks)/pow(Fum*v17_Kp + Mal*v17_Ks + v17_Kp*v17_Ks, 2);
            break;
        case 69:
            dwdp[11] = 1.0*Fum*v17_KcF*v17_v17_FM/(Fum*v17_Kp + Mal*v17_Ks + v17_Kp*v17_Ks) + 1.0*v17_v17_FM*(-Fum - v17_Ks)*(Fum*v17_KcF*v17_Kp - Mal*v17_KcR*v17_Ks)/pow(Fum*v17_Kp + Mal*v17_Ks + v17_Kp*v17_Ks, 2);
            break;
        case 70:
            dwdp[12] = 1.0*(Mal*NAD_p*v18_KcF/(v18_KiS1*v18_KmS2) - NADH*OXA*v18_KcR/(v18_KiP2*v18_KmP1))/(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*v18_KmS2) + Mal*NAD_p/(v18_KiS1*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*v18_KmP1) + Mal/v18_KiS1 + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*v18_KmS2) + NADH*OXA*v18_KiP2/v18_KmP1 + NADH/v18_KiP2 + NAD_p*v18_KmS1/(v18_KiS1*v18_KmS2) + OXA*v18_KmP2/(v18_KiP2*v18_KmP1) + 1);
            break;
        case 71:
            dwdp[12] = -1.0*NADH*OXA*v18_v18_MDH/(v18_KiP2*v18_KmP1*(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*v18_KmS2) + Mal*NAD_p/(v18_KiS1*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*v18_KmP1) + Mal/v18_KiS1 + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*v18_KmS2) + NADH*OXA*v18_KiP2/v18_KmP1 + NADH/v18_KiP2 + NAD_p*v18_KmS1/(v18_KiS1*v18_KmS2) + OXA*v18_KmP2/(v18_KiP2*v18_KmP1) + 1));
            break;
        case 72:
            dwdp[12] = 1.0*Mal*NAD_p*v18_v18_MDH/(v18_KiS1*v18_KmS2*(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*v18_KmS2) + Mal*NAD_p/(v18_KiS1*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*v18_KmP1) + Mal/v18_KiS1 + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*v18_KmS2) + NADH*OXA*v18_KiP2/v18_KmP1 + NADH/v18_KiP2 + NAD_p*v18_KmS1/(v18_KiS1*v18_KmS2) + OXA*v18_KmP2/(v18_KiP2*v18_KmP1) + 1));
            break;
        case 73:
            dwdp[12] = 1.0*NADH*OXA*v18_KcR*v18_v18_MDH/(pow(v18_KiP2, 2)*v18_KmP1*(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*v18_KmS2) + Mal*NAD_p/(v18_KiS1*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*v18_KmP1) + Mal/v18_KiS1 + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*v18_KmS2) + NADH*OXA*v18_KiP2/v18_KmP1 + NADH/v18_KiP2 + NAD_p*v18_KmS1/(v18_KiS1*v18_KmS2) + OXA*v18_KmP2/(v18_KiP2*v18_KmP1) + 1)) + 1.0*v18_v18_MDH*(Mal*NAD_p*v18_KcF/(v18_KiS1*v18_KmS2) - NADH*OXA*v18_KcR/(v18_KiP2*v18_KmP1))*(Mal*OXA*v18_KmP2/(pow(v18_KiP2, 2)*v18_KiS1*v18_KmP1) + NADH*NAD_p*OXA/(pow(v18_KiP2, 2)*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(pow(v18_KiP2, 2)*v18_KiS1*v18_KmS2) - NADH*OXA/v18_KmP1 + NADH/pow(v18_KiP2, 2) + OXA*v18_KmP2/(pow(v18_KiP2, 2)*v18_KmP1))/pow(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*v18_KmS2) + Mal*NAD_p/(v18_KiS1*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*v18_KmP1) + Mal/v18_KiS1 + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*v18_KmS2) + NADH*OXA*v18_KiP2/v18_KmP1 + NADH/v18_KiP2 + NAD_p*v18_KmS1/(v18_KiS1*v18_KmS2) + OXA*v18_KmP2/(v18_KiP2*v18_KmP1) + 1, 2);
            break;
        case 74:
            dwdp[12] = 1.0*Mal*NAD_p*OXA*v18_v18_MDH*(Mal*NAD_p*v18_KcF/(v18_KiS1*v18_KmS2) - NADH*OXA*v18_KcR/(v18_KiP2*v18_KmP1))/(pow(v18_KiP1, 2)*v18_KiS1*v18_KmS2*pow(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*v18_KmS2) + Mal*NAD_p/(v18_KiS1*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*v18_KmP1) + Mal/v18_KiS1 + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*v18_KmS2) + NADH*OXA*v18_KiP2/v18_KmP1 + NADH/v18_KiP2 + NAD_p*v18_KmS1/(v18_KiS1*v18_KmS2) + OXA*v18_KmP2/(v18_KiP2*v18_KmP1) + 1, 2));
            break;
        case 75:
            dwdp[12] = 1.0*NADH*NAD_p*OXA*v18_v18_MDH*(Mal*NAD_p*v18_KcF/(v18_KiS1*v18_KmS2) - NADH*OXA*v18_KcR/(v18_KiP2*v18_KmP1))/(v18_KiP2*pow(v18_KiS2, 2)*v18_KmP1*pow(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*v18_KmS2) + Mal*NAD_p/(v18_KiS1*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*v18_KmP1) + Mal/v18_KiS1 + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*v18_KmS2) + NADH*OXA*v18_KiP2/v18_KmP1 + NADH/v18_KiP2 + NAD_p*v18_KmS1/(v18_KiS1*v18_KmS2) + OXA*v18_KmP2/(v18_KiP2*v18_KmP1) + 1, 2));
            break;
        case 76:
            dwdp[12] = -1.0*Mal*NAD_p*v18_KcF*v18_v18_MDH/(pow(v18_KiS1, 2)*v18_KmS2*(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*v18_KmS2) + Mal*NAD_p/(v18_KiS1*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*v18_KmP1) + Mal/v18_KiS1 + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*v18_KmS2) + NADH*OXA*v18_KiP2/v18_KmP1 + NADH/v18_KiP2 + NAD_p*v18_KmS1/(v18_KiS1*v18_KmS2) + OXA*v18_KmP2/(v18_KiP2*v18_KmP1) + 1)) + 1.0*v18_v18_MDH*(Mal*NAD_p*v18_KcF/(v18_KiS1*v18_KmS2) - NADH*OXA*v18_KcR/(v18_KiP2*v18_KmP1))*(Mal*NAD_p*OXA/(v18_KiP1*pow(v18_KiS1, 2)*v18_KmS2) + Mal*NAD_p/(pow(v18_KiS1, 2)*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*pow(v18_KiS1, 2)*v18_KmP1) + Mal/pow(v18_KiS1, 2) + NADH*NAD_p*v18_KmS1/(v18_KiP2*pow(v18_KiS1, 2)*v18_KmS2) + NAD_p*v18_KmS1/(pow(v18_KiS1, 2)*v18_KmS2))/pow(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*v18_KmS2) + Mal*NAD_p/(v18_KiS1*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*v18_KmP1) + Mal/v18_KiS1 + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*v18_KmS2) + NADH*OXA*v18_KiP2/v18_KmP1 + NADH/v18_KiP2 + NAD_p*v18_KmS1/(v18_KiS1*v18_KmS2) + OXA*v18_KmP2/(v18_KiP2*v18_KmP1) + 1, 2);
            break;
        case 77:
            dwdp[12] = 1.0*v18_v18_MDH*(Mal*NAD_p*v18_KcF/(v18_KiS1*v18_KmS2) - NADH*OXA*v18_KcR/(v18_KiP2*v18_KmP1))*(-Mal*OXA/(v18_KiP2*v18_KiS1*v18_KmP1) - OXA/(v18_KiP2*v18_KmP1))/pow(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*v18_KmS2) + Mal*NAD_p/(v18_KiS1*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*v18_KmP1) + Mal/v18_KiS1 + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*v18_KmS2) + NADH*OXA*v18_KiP2/v18_KmP1 + NADH/v18_KiP2 + NAD_p*v18_KmS1/(v18_KiS1*v18_KmS2) + OXA*v18_KmP2/(v18_KiP2*v18_KmP1) + 1, 2);
            break;
        case 78:
            dwdp[12] = 1.0*NADH*OXA*v18_KcR*v18_v18_MDH/(v18_KiP2*pow(v18_KmP1, 2)*(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*v18_KmS2) + Mal*NAD_p/(v18_KiS1*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*v18_KmP1) + Mal/v18_KiS1 + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*v18_KmS2) + NADH*OXA*v18_KiP2/v18_KmP1 + NADH/v18_KiP2 + NAD_p*v18_KmS1/(v18_KiS1*v18_KmS2) + OXA*v18_KmP2/(v18_KiP2*v18_KmP1) + 1)) + 1.0*v18_v18_MDH*(Mal*NAD_p*v18_KcF/(v18_KiS1*v18_KmS2) - NADH*OXA*v18_KcR/(v18_KiP2*v18_KmP1))*(Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*pow(v18_KmP1, 2)) + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*pow(v18_KmP1, 2)) + NADH*OXA*v18_KiP2/pow(v18_KmP1, 2) + OXA*v18_KmP2/(v18_KiP2*pow(v18_KmP1, 2)))/pow(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*v18_KmS2) + Mal*NAD_p/(v18_KiS1*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*v18_KmP1) + Mal/v18_KiS1 + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*v18_KmS2) + NADH*OXA*v18_KiP2/v18_KmP1 + NADH/v18_KiP2 + NAD_p*v18_KmS1/(v18_KiS1*v18_KmS2) + OXA*v18_KmP2/(v18_KiP2*v18_KmP1) + 1, 2);
            break;
        case 79:
            dwdp[12] = -1.0*Mal*NAD_p*v18_KcF*v18_v18_MDH/(v18_KiS1*pow(v18_KmS2, 2)*(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*v18_KmS2) + Mal*NAD_p/(v18_KiS1*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*v18_KmP1) + Mal/v18_KiS1 + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*v18_KmS2) + NADH*OXA*v18_KiP2/v18_KmP1 + NADH/v18_KiP2 + NAD_p*v18_KmS1/(v18_KiS1*v18_KmS2) + OXA*v18_KmP2/(v18_KiP2*v18_KmP1) + 1)) + 1.0*v18_v18_MDH*(Mal*NAD_p*v18_KcF/(v18_KiS1*v18_KmS2) - NADH*OXA*v18_KcR/(v18_KiP2*v18_KmP1))*(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*pow(v18_KmS2, 2)) + Mal*NAD_p/(v18_KiS1*pow(v18_KmS2, 2)) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*pow(v18_KmS2, 2)) + NAD_p*v18_KmS1/(v18_KiS1*pow(v18_KmS2, 2)))/pow(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*v18_KmS2) + Mal*NAD_p/(v18_KiS1*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*v18_KmP1) + Mal/v18_KiS1 + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*v18_KmS2) + NADH*OXA*v18_KiP2/v18_KmP1 + NADH/v18_KiP2 + NAD_p*v18_KmS1/(v18_KiS1*v18_KmS2) + OXA*v18_KmP2/(v18_KiP2*v18_KmP1) + 1, 2);
            break;
        case 80:
            dwdp[12] = 1.0*v18_v18_MDH*(Mal*NAD_p*v18_KcF/(v18_KiS1*v18_KmS2) - NADH*OXA*v18_KcR/(v18_KiP2*v18_KmP1))*(-NADH*NAD_p/(v18_KiP2*v18_KiS1*v18_KmS2) - NAD_p/(v18_KiS1*v18_KmS2))/pow(Mal*NAD_p*OXA/(v18_KiP1*v18_KiS1*v18_KmS2) + Mal*NAD_p/(v18_KiS1*v18_KmS2) + Mal*OXA*v18_KmP2/(v18_KiP2*v18_KiS1*v18_KmP1) + Mal/v18_KiS1 + NADH*NAD_p*OXA/(v18_KiP2*v18_KiS2*v18_KmP1) + NADH*NAD_p*v18_KmS1/(v18_KiP2*v18_KiS1*v18_KmS2) + NADH*OXA*v18_KiP2/v18_KmP1 + NADH/v18_KiP2 + NAD_p*v18_KmS1/(v18_KiS1*v18_KmS2) + OXA*v18_KmP2/(v18_KiP2*v18_KmP1) + 1, 2);
            break;
        case 81:
            dwdp[13] = -1.0*ATP_cyt*pow(F6P, 2)*v2_V2/(pow(ATP_cyt + v2_K2ATP, 2)*(pow(F6P, 2) + v2_K2*(1 + pow(ATP_cyt, 2)*v2_k2/pow(AMP, 2))));
            break;
        case 82:
            dwdp[13] = -1.0*pow(ATP_cyt, 3)*pow(F6P, 2)*v2_K2*v2_V2/(pow(AMP, 2)*(ATP_cyt + v2_K2ATP)*pow(pow(F6P, 2) + v2_K2*(1 + pow(ATP_cyt, 2)*v2_k2/pow(AMP, 2)), 2));
            break;
        case 83:
            dwdp[13] = 1.0*ATP_cyt*pow(F6P, 2)*v2_V2*(-1 - pow(ATP_cyt, 2)*v2_k2/pow(AMP, 2))/((ATP_cyt + v2_K2ATP)*pow(pow(F6P, 2) + v2_K2*(1 + pow(ATP_cyt, 2)*v2_k2/pow(AMP, 2)), 2));
            break;
        case 84:
            dwdp[13] = 1.0*ATP_cyt*pow(F6P, 2)/((ATP_cyt + v2_K2ATP)*(pow(F6P, 2) + v2_K2*(1 + pow(ATP_cyt, 2)*v2_k2/pow(AMP, 2))));
            break;
        case 85:
            dwdp[14] = 1.0*v20_KcF*v20_KcR*(Ala*OG - Glu*Pyr/v20_Keq)/(Ala*Glu*v20_KcF*v20_KmP2/(v20_Keq*v20_KiS1) + Ala*OG*v20_KcR + Ala*v20_KcR*v20_KmS2 + Glu*Pyr*v20_KcF/v20_Keq + Glu*v20_KcF*v20_KmP2/v20_Keq + OG*Pyr*v20_KcR*v20_KmS1/v20_KiP2 + OG*v20_KcR*v20_KmS1 + Pyr*v20_KcF*v20_KmP1/v20_Keq);
            break;
        case 86:
            dwdp[14] = 1.0*v20_KcF*v20_KcR*v20_v20_AlaTA*(Ala*OG - Glu*Pyr/v20_Keq)*(-Ala*OG - Ala*v20_KmS2 - OG*Pyr*v20_KmS1/v20_KiP2 - OG*v20_KmS1)/pow(Ala*Glu*v20_KcF*v20_KmP2/(v20_Keq*v20_KiS1) + Ala*OG*v20_KcR + Ala*v20_KcR*v20_KmS2 + Glu*Pyr*v20_KcF/v20_Keq + Glu*v20_KcF*v20_KmP2/v20_Keq + OG*Pyr*v20_KcR*v20_KmS1/v20_KiP2 + OG*v20_KcR*v20_KmS1 + Pyr*v20_KcF*v20_KmP1/v20_Keq, 2) + 1.0*v20_KcF*v20_v20_AlaTA*(Ala*OG - Glu*Pyr/v20_Keq)/(Ala*Glu*v20_KcF*v20_KmP2/(v20_Keq*v20_KiS1) + Ala*OG*v20_KcR + Ala*v20_KcR*v20_KmS2 + Glu*Pyr*v20_KcF/v20_Keq + Glu*v20_KcF*v20_KmP2/v20_Keq + OG*Pyr*v20_KcR*v20_KmS1/v20_KiP2 + OG*v20_KcR*v20_KmS1 + Pyr*v20_KcF*v20_KmP1/v20_Keq);
            break;
        case 87:
            dwdp[14] = 1.0*v20_KcF*v20_KcR*v20_v20_AlaTA*(Ala*OG - Glu*Pyr/v20_Keq)*(-Ala*Glu*v20_KmP2/(v20_Keq*v20_KiS1) - Glu*Pyr/v20_Keq - Glu*v20_KmP2/v20_Keq - Pyr*v20_KmP1/v20_Keq)/pow(Ala*Glu*v20_KcF*v20_KmP2/(v20_Keq*v20_KiS1) + Ala*OG*v20_KcR + Ala*v20_KcR*v20_KmS2 + Glu*Pyr*v20_KcF/v20_Keq + Glu*v20_KcF*v20_KmP2/v20_Keq + OG*Pyr*v20_KcR*v20_KmS1/v20_KiP2 + OG*v20_KcR*v20_KmS1 + Pyr*v20_KcF*v20_KmP1/v20_Keq, 2) + 1.0*v20_KcR*v20_v20_AlaTA*(Ala*OG - Glu*Pyr/v20_Keq)/(Ala*Glu*v20_KcF*v20_KmP2/(v20_Keq*v20_KiS1) + Ala*OG*v20_KcR + Ala*v20_KcR*v20_KmS2 + Glu*Pyr*v20_KcF/v20_Keq + Glu*v20_KcF*v20_KmP2/v20_Keq + OG*Pyr*v20_KcR*v20_KmS1/v20_KiP2 + OG*v20_KcR*v20_KmS1 + Pyr*v20_KcF*v20_KmP1/v20_Keq);
            break;
        case 88:
            dwdp[14] = 1.0*Glu*Pyr*v20_KcF*v20_KcR*v20_v20_AlaTA/(pow(v20_Keq, 2)*(Ala*Glu*v20_KcF*v20_KmP2/(v20_Keq*v20_KiS1) + Ala*OG*v20_KcR + Ala*v20_KcR*v20_KmS2 + Glu*Pyr*v20_KcF/v20_Keq + Glu*v20_KcF*v20_KmP2/v20_Keq + OG*Pyr*v20_KcR*v20_KmS1/v20_KiP2 + OG*v20_KcR*v20_KmS1 + Pyr*v20_KcF*v20_KmP1/v20_Keq)) + 1.0*v20_KcF*v20_KcR*v20_v20_AlaTA*(Ala*OG - Glu*Pyr/v20_Keq)*(Ala*Glu*v20_KcF*v20_KmP2/(pow(v20_Keq, 2)*v20_KiS1) + Glu*Pyr*v20_KcF/pow(v20_Keq, 2) + Glu*v20_KcF*v20_KmP2/pow(v20_Keq, 2) + Pyr*v20_KcF*v20_KmP1/pow(v20_Keq, 2))/pow(Ala*Glu*v20_KcF*v20_KmP2/(v20_Keq*v20_KiS1) + Ala*OG*v20_KcR + Ala*v20_KcR*v20_KmS2 + Glu*Pyr*v20_KcF/v20_Keq + Glu*v20_KcF*v20_KmP2/v20_Keq + OG*Pyr*v20_KcR*v20_KmS1/v20_KiP2 + OG*v20_KcR*v20_KmS1 + Pyr*v20_KcF*v20_KmP1/v20_Keq, 2);
            break;
        case 89:
            dwdp[14] = 1.0*OG*Pyr*v20_KcF*pow(v20_KcR, 2)*v20_KmS1*v20_v20_AlaTA*(Ala*OG - Glu*Pyr/v20_Keq)/(pow(v20_KiP2, 2)*pow(Ala*Glu*v20_KcF*v20_KmP2/(v20_Keq*v20_KiS1) + Ala*OG*v20_KcR + Ala*v20_KcR*v20_KmS2 + Glu*Pyr*v20_KcF/v20_Keq + Glu*v20_KcF*v20_KmP2/v20_Keq + OG*Pyr*v20_KcR*v20_KmS1/v20_KiP2 + OG*v20_KcR*v20_KmS1 + Pyr*v20_KcF*v20_KmP1/v20_Keq, 2));
            break;
        case 90:
            dwdp[14] = 1.0*Ala*Glu*pow(v20_KcF, 2)*v20_KcR*v20_KmP2*v20_v20_AlaTA*(Ala*OG - Glu*Pyr/v20_Keq)/(v20_Keq*pow(v20_KiS1, 2)*pow(Ala*Glu*v20_KcF*v20_KmP2/(v20_Keq*v20_KiS1) + Ala*OG*v20_KcR + Ala*v20_KcR*v20_KmS2 + Glu*Pyr*v20_KcF/v20_Keq + Glu*v20_KcF*v20_KmP2/v20_Keq + OG*Pyr*v20_KcR*v20_KmS1/v20_KiP2 + OG*v20_KcR*v20_KmS1 + Pyr*v20_KcF*v20_KmP1/v20_Keq, 2));
            break;
        case 91:
            dwdp[14] = 1.0*v20_KcF*v20_KcR*v20_v20_AlaTA*(Ala*OG - Glu*Pyr/v20_Keq)*(-Ala*Glu*v20_KcF/(v20_Keq*v20_KiS1) - Glu*v20_KcF/v20_Keq)/pow(Ala*Glu*v20_KcF*v20_KmP2/(v20_Keq*v20_KiS1) + Ala*OG*v20_KcR + Ala*v20_KcR*v20_KmS2 + Glu*Pyr*v20_KcF/v20_Keq + Glu*v20_KcF*v20_KmP2/v20_Keq + OG*Pyr*v20_KcR*v20_KmS1/v20_KiP2 + OG*v20_KcR*v20_KmS1 + Pyr*v20_KcF*v20_KmP1/v20_Keq, 2);
            break;
        case 92:
            dwdp[14] = -1.0*Pyr*pow(v20_KcF, 2)*v20_KcR*v20_v20_AlaTA*(Ala*OG - Glu*Pyr/v20_Keq)/(v20_Keq*pow(Ala*Glu*v20_KcF*v20_KmP2/(v20_Keq*v20_KiS1) + Ala*OG*v20_KcR + Ala*v20_KcR*v20_KmS2 + Glu*Pyr*v20_KcF/v20_Keq + Glu*v20_KcF*v20_KmP2/v20_Keq + OG*Pyr*v20_KcR*v20_KmS1/v20_KiP2 + OG*v20_KcR*v20_KmS1 + Pyr*v20_KcF*v20_KmP1/v20_Keq, 2));
            break;
        case 93:
            dwdp[14] = -1.0*Ala*v20_KcF*pow(v20_KcR, 2)*v20_v20_AlaTA*(Ala*OG - Glu*Pyr/v20_Keq)/pow(Ala*Glu*v20_KcF*v20_KmP2/(v20_Keq*v20_KiS1) + Ala*OG*v20_KcR + Ala*v20_KcR*v20_KmS2 + Glu*Pyr*v20_KcF/v20_Keq + Glu*v20_KcF*v20_KmP2/v20_Keq + OG*Pyr*v20_KcR*v20_KmS1/v20_KiP2 + OG*v20_KcR*v20_KmS1 + Pyr*v20_KcF*v20_KmP1/v20_Keq, 2);
            break;
        case 94:
            dwdp[14] = 1.0*v20_KcF*v20_KcR*v20_v20_AlaTA*(Ala*OG - Glu*Pyr/v20_Keq)*(-OG*Pyr*v20_KcR/v20_KiP2 - OG*v20_KcR)/pow(Ala*Glu*v20_KcF*v20_KmP2/(v20_Keq*v20_KiS1) + Ala*OG*v20_KcR + Ala*v20_KcR*v20_KmS2 + Glu*Pyr*v20_KcF/v20_Keq + Glu*v20_KcF*v20_KmP2/v20_Keq + OG*Pyr*v20_KcR*v20_KmS1/v20_KiP2 + OG*v20_KcR*v20_KmS1 + Pyr*v20_KcF*v20_KmP1/v20_Keq, 2);
            break;
        case 95:
            dwdp[15] = 1.0*v21_KcF*v21_KcR*(-Asp*OG/v21_Keq + Glu*OXA)/(Asp*OG*v21_KcF/v21_Keq + Asp*OXA*v21_KcF*v21_KmP2/(v21_Keq*v21_KiS1) + Asp*v21_KcF*v21_KmP2/v21_Keq + Glu*OG*v21_KcR*v21_KmS1/v21_KiP2 + Glu*OXA*v21_KcR + Glu*v21_KcR*v21_KmS1 + OG*v21_KcF*v21_KmP1/v21_Keq + OXA*v21_KcR*v21_KmS2);
            break;
        case 96:
            dwdp[15] = 1.0*v21_KcF*v21_KcR*v21_v21_AspTA*(-Asp*OG/v21_Keq + Glu*OXA)*(-Glu*OG*v21_KmS1/v21_KiP2 - Glu*OXA - Glu*v21_KmS1 - OXA*v21_KmS2)/pow(Asp*OG*v21_KcF/v21_Keq + Asp*OXA*v21_KcF*v21_KmP2/(v21_Keq*v21_KiS1) + Asp*v21_KcF*v21_KmP2/v21_Keq + Glu*OG*v21_KcR*v21_KmS1/v21_KiP2 + Glu*OXA*v21_KcR + Glu*v21_KcR*v21_KmS1 + OG*v21_KcF*v21_KmP1/v21_Keq + OXA*v21_KcR*v21_KmS2, 2) + 1.0*v21_KcF*v21_v21_AspTA*(-Asp*OG/v21_Keq + Glu*OXA)/(Asp*OG*v21_KcF/v21_Keq + Asp*OXA*v21_KcF*v21_KmP2/(v21_Keq*v21_KiS1) + Asp*v21_KcF*v21_KmP2/v21_Keq + Glu*OG*v21_KcR*v21_KmS1/v21_KiP2 + Glu*OXA*v21_KcR + Glu*v21_KcR*v21_KmS1 + OG*v21_KcF*v21_KmP1/v21_Keq + OXA*v21_KcR*v21_KmS2);
            break;
        case 97:
            dwdp[15] = 1.0*v21_KcF*v21_KcR*v21_v21_AspTA*(-Asp*OG/v21_Keq + Glu*OXA)*(-Asp*OG/v21_Keq - Asp*OXA*v21_KmP2/(v21_Keq*v21_KiS1) - Asp*v21_KmP2/v21_Keq - OG*v21_KmP1/v21_Keq)/pow(Asp*OG*v21_KcF/v21_Keq + Asp*OXA*v21_KcF*v21_KmP2/(v21_Keq*v21_KiS1) + Asp*v21_KcF*v21_KmP2/v21_Keq + Glu*OG*v21_KcR*v21_KmS1/v21_KiP2 + Glu*OXA*v21_KcR + Glu*v21_KcR*v21_KmS1 + OG*v21_KcF*v21_KmP1/v21_Keq + OXA*v21_KcR*v21_KmS2, 2) + 1.0*v21_KcR*v21_v21_AspTA*(-Asp*OG/v21_Keq + Glu*OXA)/(Asp*OG*v21_KcF/v21_Keq + Asp*OXA*v21_KcF*v21_KmP2/(v21_Keq*v21_KiS1) + Asp*v21_KcF*v21_KmP2/v21_Keq + Glu*OG*v21_KcR*v21_KmS1/v21_KiP2 + Glu*OXA*v21_KcR + Glu*v21_KcR*v21_KmS1 + OG*v21_KcF*v21_KmP1/v21_Keq + OXA*v21_KcR*v21_KmS2);
            break;
        case 98:
            dwdp[15] = 1.0*Asp*OG*v21_KcF*v21_KcR*v21_v21_AspTA/(pow(v21_Keq, 2)*(Asp*OG*v21_KcF/v21_Keq + Asp*OXA*v21_KcF*v21_KmP2/(v21_Keq*v21_KiS1) + Asp*v21_KcF*v21_KmP2/v21_Keq + Glu*OG*v21_KcR*v21_KmS1/v21_KiP2 + Glu*OXA*v21_KcR + Glu*v21_KcR*v21_KmS1 + OG*v21_KcF*v21_KmP1/v21_Keq + OXA*v21_KcR*v21_KmS2)) + 1.0*v21_KcF*v21_KcR*v21_v21_AspTA*(-Asp*OG/v21_Keq + Glu*OXA)*(Asp*OG*v21_KcF/pow(v21_Keq, 2) + Asp*OXA*v21_KcF*v21_KmP2/(pow(v21_Keq, 2)*v21_KiS1) + Asp*v21_KcF*v21_KmP2/pow(v21_Keq, 2) + OG*v21_KcF*v21_KmP1/pow(v21_Keq, 2))/pow(Asp*OG*v21_KcF/v21_Keq + Asp*OXA*v21_KcF*v21_KmP2/(v21_Keq*v21_KiS1) + Asp*v21_KcF*v21_KmP2/v21_Keq + Glu*OG*v21_KcR*v21_KmS1/v21_KiP2 + Glu*OXA*v21_KcR + Glu*v21_KcR*v21_KmS1 + OG*v21_KcF*v21_KmP1/v21_Keq + OXA*v21_KcR*v21_KmS2, 2);
            break;
        case 99:
            dwdp[15] = 1.0*Glu*OG*v21_KcF*pow(v21_KcR, 2)*v21_KmS1*v21_v21_AspTA*(-Asp*OG/v21_Keq + Glu*OXA)/(pow(v21_KiP2, 2)*pow(Asp*OG*v21_KcF/v21_Keq + Asp*OXA*v21_KcF*v21_KmP2/(v21_Keq*v21_KiS1) + Asp*v21_KcF*v21_KmP2/v21_Keq + Glu*OG*v21_KcR*v21_KmS1/v21_KiP2 + Glu*OXA*v21_KcR + Glu*v21_KcR*v21_KmS1 + OG*v21_KcF*v21_KmP1/v21_Keq + OXA*v21_KcR*v21_KmS2, 2));
            break;
        case 100:
            dwdp[15] = 1.0*Asp*OXA*pow(v21_KcF, 2)*v21_KcR*v21_KmP2*v21_v21_AspTA*(-Asp*OG/v21_Keq + Glu*OXA)/(v21_Keq*pow(v21_KiS1, 2)*pow(Asp*OG*v21_KcF/v21_Keq + Asp*OXA*v21_KcF*v21_KmP2/(v21_Keq*v21_KiS1) + Asp*v21_KcF*v21_KmP2/v21_Keq + Glu*OG*v21_KcR*v21_KmS1/v21_KiP2 + Glu*OXA*v21_KcR + Glu*v21_KcR*v21_KmS1 + OG*v21_KcF*v21_KmP1/v21_Keq + OXA*v21_KcR*v21_KmS2, 2));
            break;
        case 101:
            dwdp[15] = 1.0*v21_KcF*v21_KcR*v21_v21_AspTA*(-Asp*OG/v21_Keq + Glu*OXA)*(-Asp*OXA*v21_KcF/(v21_Keq*v21_KiS1) - Asp*v21_KcF/v21_Keq)/pow(Asp*OG*v21_KcF/v21_Keq + Asp*OXA*v21_KcF*v21_KmP2/(v21_Keq*v21_KiS1) + Asp*v21_KcF*v21_KmP2/v21_Keq + Glu*OG*v21_KcR*v21_KmS1/v21_KiP2 + Glu*OXA*v21_KcR + Glu*v21_KcR*v21_KmS1 + OG*v21_KcF*v21_KmP1/v21_Keq + OXA*v21_KcR*v21_KmS2, 2);
            break;
        case 102:
            dwdp[15] = -1.0*OG*pow(v21_KcF, 2)*v21_KcR*v21_v21_AspTA*(-Asp*OG/v21_Keq + Glu*OXA)/(v21_Keq*pow(Asp*OG*v21_KcF/v21_Keq + Asp*OXA*v21_KcF*v21_KmP2/(v21_Keq*v21_KiS1) + Asp*v21_KcF*v21_KmP2/v21_Keq + Glu*OG*v21_KcR*v21_KmS1/v21_KiP2 + Glu*OXA*v21_KcR + Glu*v21_KcR*v21_KmS1 + OG*v21_KcF*v21_KmP1/v21_Keq + OXA*v21_KcR*v21_KmS2, 2));
            break;
        case 103:
            dwdp[15] = -1.0*OXA*v21_KcF*pow(v21_KcR, 2)*v21_v21_AspTA*(-Asp*OG/v21_Keq + Glu*OXA)/pow(Asp*OG*v21_KcF/v21_Keq + Asp*OXA*v21_KcF*v21_KmP2/(v21_Keq*v21_KiS1) + Asp*v21_KcF*v21_KmP2/v21_Keq + Glu*OG*v21_KcR*v21_KmS1/v21_KiP2 + Glu*OXA*v21_KcR + Glu*v21_KcR*v21_KmS1 + OG*v21_KcF*v21_KmP1/v21_Keq + OXA*v21_KcR*v21_KmS2, 2);
            break;
        case 104:
            dwdp[15] = 1.0*v21_KcF*v21_KcR*v21_v21_AspTA*(-Asp*OG/v21_Keq + Glu*OXA)*(-Glu*OG*v21_KcR/v21_KiP2 - Glu*v21_KcR)/pow(Asp*OG*v21_KcF/v21_Keq + Asp*OXA*v21_KcF*v21_KmP2/(v21_Keq*v21_KiS1) + Asp*v21_KcF*v21_KmP2/v21_Keq + Glu*OG*v21_KcR*v21_KmS1/v21_KiP2 + Glu*OXA*v21_KcR + Glu*v21_KcR*v21_KmS1 + OG*v21_KcF*v21_KmP1/v21_Keq + OXA*v21_KcR*v21_KmS2, 2);
            break;
        case 105:
            dwdp[16] = 1.0*(Asp*Glu_cyt*v22_KcF/(v22_KiS1*v22_KiS2*v22_alpha) - Asp_cyt*Glu*v22_KcR/(v22_KiP1*v22_KiP2*v22_beta))/(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1);
            break;
        case 106:
            dwdp[16] = 1.0*Asp*Glu*v22_v22_AGC*(Asp*Glu_cyt*v22_KcF/(v22_KiS1*v22_KiS2*v22_alpha) - Asp_cyt*Glu*v22_KcR/(v22_KiP1*v22_KiP2*v22_beta))/(v22_KiP1*v22_KiS1*pow(v22_delta, 2)*pow(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1, 2));
            break;
        case 107:
            dwdp[16] = 1.0*Asp_cyt*Glu_cyt*v22_v22_AGC*(Asp*Glu_cyt*v22_KcF/(v22_KiS1*v22_KiS2*v22_alpha) - Asp_cyt*Glu*v22_KcR/(v22_KiP1*v22_KiP2*v22_beta))/(v22_KiP2*v22_KiS2*pow(v22_gamma, 2)*pow(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1, 2));
            break;
        case 108:
            dwdp[16] = 1.0*Asp_cyt*Glu*v22_KcR*v22_v22_AGC/(v22_KiP1*v22_KiP2*pow(v22_beta, 2)*(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1)) + 1.0*Asp_cyt*Glu*v22_v22_AGC*(Asp*Glu_cyt*v22_KcF/(v22_KiS1*v22_KiS2*v22_alpha) - Asp_cyt*Glu*v22_KcR/(v22_KiP1*v22_KiP2*v22_beta))/(v22_KiP1*v22_KiP2*pow(v22_beta, 2)*pow(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1, 2));
            break;
        case 109:
            dwdp[16] = -1.0*Asp*Glu_cyt*v22_KcF*v22_v22_AGC/(v22_KiS1*v22_KiS2*pow(v22_alpha, 2)*(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1)) + 1.0*Asp*Glu_cyt*v22_v22_AGC*(Asp*Glu_cyt*v22_KcF/(v22_KiS1*v22_KiS2*v22_alpha) - Asp_cyt*Glu*v22_KcR/(v22_KiP1*v22_KiP2*v22_beta))/(v22_KiS1*v22_KiS2*pow(v22_alpha, 2)*pow(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1, 2));
            break;
        case 110:
            dwdp[16] = -1.0*Asp_cyt*Glu*v22_v22_AGC/(v22_KiP1*v22_KiP2*v22_beta*(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1));
            break;
        case 111:
            dwdp[16] = 1.0*Asp*Glu_cyt*v22_v22_AGC/(v22_KiS1*v22_KiS2*v22_alpha*(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1));
            break;
        case 112:
            dwdp[16] = 1.0*Asp_cyt*Glu*v22_KcR*v22_v22_AGC/(v22_KiP1*pow(v22_KiP2, 2)*v22_beta*(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1)) + 1.0*v22_v22_AGC*(Asp*Glu_cyt*v22_KcF/(v22_KiS1*v22_KiS2*v22_alpha) - Asp_cyt*Glu*v22_KcR/(v22_KiP1*v22_KiP2*v22_beta))*(Asp_cyt*Glu/(v22_KiP1*pow(v22_KiP2, 2)*v22_beta) + Asp_cyt*Glu_cyt/(pow(v22_KiP2, 2)*v22_KiS2*v22_gamma) + Asp_cyt/pow(v22_KiP2, 2))/pow(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1, 2);
            break;
        case 113:
            dwdp[16] = 1.0*Asp_cyt*Glu*v22_KcR*v22_v22_AGC/(pow(v22_KiP1, 2)*v22_KiP2*v22_beta*(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1)) + 1.0*v22_v22_AGC*(Asp*Glu_cyt*v22_KcF/(v22_KiS1*v22_KiS2*v22_alpha) - Asp_cyt*Glu*v22_KcR/(v22_KiP1*v22_KiP2*v22_beta))*(Asp*Glu/(pow(v22_KiP1, 2)*v22_KiS1*v22_delta) + Asp_cyt*Glu/(pow(v22_KiP1, 2)*v22_KiP2*v22_beta) + Glu/pow(v22_KiP1, 2))/pow(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1, 2);
            break;
        case 114:
            dwdp[16] = -1.0*Asp*Glu_cyt*v22_KcF*v22_v22_AGC/(v22_KiS1*pow(v22_KiS2, 2)*v22_alpha*(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1)) + 1.0*v22_v22_AGC*(Asp*Glu_cyt*v22_KcF/(v22_KiS1*v22_KiS2*v22_alpha) - Asp_cyt*Glu*v22_KcR/(v22_KiP1*v22_KiP2*v22_beta))*(Asp*Glu_cyt/(v22_KiS1*pow(v22_KiS2, 2)*v22_alpha) + Asp_cyt*Glu_cyt/(v22_KiP2*pow(v22_KiS2, 2)*v22_gamma) + Glu_cyt/pow(v22_KiS2, 2))/pow(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1, 2);
            break;
        case 115:
            dwdp[16] = -1.0*Asp*Glu_cyt*v22_KcF*v22_v22_AGC/(pow(v22_KiS1, 2)*v22_KiS2*v22_alpha*(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1)) + 1.0*v22_v22_AGC*(Asp*Glu_cyt*v22_KcF/(v22_KiS1*v22_KiS2*v22_alpha) - Asp_cyt*Glu*v22_KcR/(v22_KiP1*v22_KiP2*v22_beta))*(Asp*Glu/(v22_KiP1*pow(v22_KiS1, 2)*v22_delta) + Asp*Glu_cyt/(pow(v22_KiS1, 2)*v22_KiS2*v22_alpha) + Asp/pow(v22_KiS1, 2))/pow(Asp*Glu/(v22_KiP1*v22_KiS1*v22_delta) + Asp*Glu_cyt/(v22_KiS1*v22_KiS2*v22_alpha) + Asp/v22_KiS1 + Asp_cyt*Glu/(v22_KiP1*v22_KiP2*v22_beta) + Asp_cyt*Glu_cyt/(v22_KiP2*v22_KiS2*v22_gamma) + Asp_cyt/v22_KiP2 + Glu/v22_KiP1 + Glu_cyt/v22_KiS2 + 1, 2);
            break;
        case 116:
            dwdp[17] = 1.0*v24_KcF*v24_KcR*(NADH*Q - NAD_p*QH2/v24_Keq)/(NADH*NAD_p*v24_KcF*v24_KmP2/(v24_Keq*v24_KiS1) + NADH*Q*v24_KcR + NADH*v24_KcR*v24_KmS2 + NAD_p*QH2*v24_KcF/v24_Keq + NAD_p*v24_KcF*v24_KmP2/v24_Keq + Q*QH2*v24_KcR*v24_KmS1/v24_KiP2 + Q*v24_KcR*v24_KmS1 + QH2*v24_KcF*v24_KmP1/v24_Keq);
            break;
        case 117:
            dwdp[17] = 1.0*v24_KcF*v24_KcR*v24_v24_Complex_I*(NADH*Q - NAD_p*QH2/v24_Keq)*(-NADH*Q - NADH*v24_KmS2 - Q*QH2*v24_KmS1/v24_KiP2 - Q*v24_KmS1)/pow(NADH*NAD_p*v24_KcF*v24_KmP2/(v24_Keq*v24_KiS1) + NADH*Q*v24_KcR + NADH*v24_KcR*v24_KmS2 + NAD_p*QH2*v24_KcF/v24_Keq + NAD_p*v24_KcF*v24_KmP2/v24_Keq + Q*QH2*v24_KcR*v24_KmS1/v24_KiP2 + Q*v24_KcR*v24_KmS1 + QH2*v24_KcF*v24_KmP1/v24_Keq, 2) + 1.0*v24_KcF*v24_v24_Complex_I*(NADH*Q - NAD_p*QH2/v24_Keq)/(NADH*NAD_p*v24_KcF*v24_KmP2/(v24_Keq*v24_KiS1) + NADH*Q*v24_KcR + NADH*v24_KcR*v24_KmS2 + NAD_p*QH2*v24_KcF/v24_Keq + NAD_p*v24_KcF*v24_KmP2/v24_Keq + Q*QH2*v24_KcR*v24_KmS1/v24_KiP2 + Q*v24_KcR*v24_KmS1 + QH2*v24_KcF*v24_KmP1/v24_Keq);
            break;
        case 118:
            dwdp[17] = 1.0*v24_KcF*v24_KcR*v24_v24_Complex_I*(NADH*Q - NAD_p*QH2/v24_Keq)*(-NADH*NAD_p*v24_KmP2/(v24_Keq*v24_KiS1) - NAD_p*QH2/v24_Keq - NAD_p*v24_KmP2/v24_Keq - QH2*v24_KmP1/v24_Keq)/pow(NADH*NAD_p*v24_KcF*v24_KmP2/(v24_Keq*v24_KiS1) + NADH*Q*v24_KcR + NADH*v24_KcR*v24_KmS2 + NAD_p*QH2*v24_KcF/v24_Keq + NAD_p*v24_KcF*v24_KmP2/v24_Keq + Q*QH2*v24_KcR*v24_KmS1/v24_KiP2 + Q*v24_KcR*v24_KmS1 + QH2*v24_KcF*v24_KmP1/v24_Keq, 2) + 1.0*v24_KcR*v24_v24_Complex_I*(NADH*Q - NAD_p*QH2/v24_Keq)/(NADH*NAD_p*v24_KcF*v24_KmP2/(v24_Keq*v24_KiS1) + NADH*Q*v24_KcR + NADH*v24_KcR*v24_KmS2 + NAD_p*QH2*v24_KcF/v24_Keq + NAD_p*v24_KcF*v24_KmP2/v24_Keq + Q*QH2*v24_KcR*v24_KmS1/v24_KiP2 + Q*v24_KcR*v24_KmS1 + QH2*v24_KcF*v24_KmP1/v24_Keq);
            break;
        case 119:
            dwdp[17] = 1.0*NAD_p*QH2*v24_KcF*v24_KcR*v24_v24_Complex_I/(pow(v24_Keq, 2)*(NADH*NAD_p*v24_KcF*v24_KmP2/(v24_Keq*v24_KiS1) + NADH*Q*v24_KcR + NADH*v24_KcR*v24_KmS2 + NAD_p*QH2*v24_KcF/v24_Keq + NAD_p*v24_KcF*v24_KmP2/v24_Keq + Q*QH2*v24_KcR*v24_KmS1/v24_KiP2 + Q*v24_KcR*v24_KmS1 + QH2*v24_KcF*v24_KmP1/v24_Keq)) + 1.0*v24_KcF*v24_KcR*v24_v24_Complex_I*(NADH*Q - NAD_p*QH2/v24_Keq)*(NADH*NAD_p*v24_KcF*v24_KmP2/(pow(v24_Keq, 2)*v24_KiS1) + NAD_p*QH2*v24_KcF/pow(v24_Keq, 2) + NAD_p*v24_KcF*v24_KmP2/pow(v24_Keq, 2) + QH2*v24_KcF*v24_KmP1/pow(v24_Keq, 2))/pow(NADH*NAD_p*v24_KcF*v24_KmP2/(v24_Keq*v24_KiS1) + NADH*Q*v24_KcR + NADH*v24_KcR*v24_KmS2 + NAD_p*QH2*v24_KcF/v24_Keq + NAD_p*v24_KcF*v24_KmP2/v24_Keq + Q*QH2*v24_KcR*v24_KmS1/v24_KiP2 + Q*v24_KcR*v24_KmS1 + QH2*v24_KcF*v24_KmP1/v24_Keq, 2);
            break;
        case 120:
            dwdp[17] = 1.0*Q*QH2*v24_KcF*pow(v24_KcR, 2)*v24_KmS1*v24_v24_Complex_I*(NADH*Q - NAD_p*QH2/v24_Keq)/(pow(v24_KiP2, 2)*pow(NADH*NAD_p*v24_KcF*v24_KmP2/(v24_Keq*v24_KiS1) + NADH*Q*v24_KcR + NADH*v24_KcR*v24_KmS2 + NAD_p*QH2*v24_KcF/v24_Keq + NAD_p*v24_KcF*v24_KmP2/v24_Keq + Q*QH2*v24_KcR*v24_KmS1/v24_KiP2 + Q*v24_KcR*v24_KmS1 + QH2*v24_KcF*v24_KmP1/v24_Keq, 2));
            break;
        case 121:
            dwdp[17] = 1.0*NADH*NAD_p*pow(v24_KcF, 2)*v24_KcR*v24_KmP2*v24_v24_Complex_I*(NADH*Q - NAD_p*QH2/v24_Keq)/(v24_Keq*pow(v24_KiS1, 2)*pow(NADH*NAD_p*v24_KcF*v24_KmP2/(v24_Keq*v24_KiS1) + NADH*Q*v24_KcR + NADH*v24_KcR*v24_KmS2 + NAD_p*QH2*v24_KcF/v24_Keq + NAD_p*v24_KcF*v24_KmP2/v24_Keq + Q*QH2*v24_KcR*v24_KmS1/v24_KiP2 + Q*v24_KcR*v24_KmS1 + QH2*v24_KcF*v24_KmP1/v24_Keq, 2));
            break;
        case 122:
            dwdp[17] = 1.0*v24_KcF*v24_KcR*v24_v24_Complex_I*(NADH*Q - NAD_p*QH2/v24_Keq)*(-NADH*NAD_p*v24_KcF/(v24_Keq*v24_KiS1) - NAD_p*v24_KcF/v24_Keq)/pow(NADH*NAD_p*v24_KcF*v24_KmP2/(v24_Keq*v24_KiS1) + NADH*Q*v24_KcR + NADH*v24_KcR*v24_KmS2 + NAD_p*QH2*v24_KcF/v24_Keq + NAD_p*v24_KcF*v24_KmP2/v24_Keq + Q*QH2*v24_KcR*v24_KmS1/v24_KiP2 + Q*v24_KcR*v24_KmS1 + QH2*v24_KcF*v24_KmP1/v24_Keq, 2);
            break;
        case 123:
            dwdp[17] = -1.0*QH2*pow(v24_KcF, 2)*v24_KcR*v24_v24_Complex_I*(NADH*Q - NAD_p*QH2/v24_Keq)/(v24_Keq*pow(NADH*NAD_p*v24_KcF*v24_KmP2/(v24_Keq*v24_KiS1) + NADH*Q*v24_KcR + NADH*v24_KcR*v24_KmS2 + NAD_p*QH2*v24_KcF/v24_Keq + NAD_p*v24_KcF*v24_KmP2/v24_Keq + Q*QH2*v24_KcR*v24_KmS1/v24_KiP2 + Q*v24_KcR*v24_KmS1 + QH2*v24_KcF*v24_KmP1/v24_Keq, 2));
            break;
        case 124:
            dwdp[17] = -1.0*NADH*v24_KcF*pow(v24_KcR, 2)*v24_v24_Complex_I*(NADH*Q - NAD_p*QH2/v24_Keq)/pow(NADH*NAD_p*v24_KcF*v24_KmP2/(v24_Keq*v24_KiS1) + NADH*Q*v24_KcR + NADH*v24_KcR*v24_KmS2 + NAD_p*QH2*v24_KcF/v24_Keq + NAD_p*v24_KcF*v24_KmP2/v24_Keq + Q*QH2*v24_KcR*v24_KmS1/v24_KiP2 + Q*v24_KcR*v24_KmS1 + QH2*v24_KcF*v24_KmP1/v24_Keq, 2);
            break;
        case 125:
            dwdp[17] = 1.0*v24_KcF*v24_KcR*v24_v24_Complex_I*(NADH*Q - NAD_p*QH2/v24_Keq)*(-Q*QH2*v24_KcR/v24_KiP2 - Q*v24_KcR)/pow(NADH*NAD_p*v24_KcF*v24_KmP2/(v24_Keq*v24_KiS1) + NADH*Q*v24_KcR + NADH*v24_KcR*v24_KmS2 + NAD_p*QH2*v24_KcF/v24_Keq + NAD_p*v24_KcF*v24_KmP2/v24_Keq + Q*QH2*v24_KcR*v24_KmS1/v24_KiP2 + Q*v24_KcR*v24_KmS1 + QH2*v24_KcF*v24_KmP1/v24_Keq, 2);
            break;
        case 126:
            dwdp[18] = 1.0*Cytc3p*QH2*v25_KcF/(Cytc2p*(Cytc3p*QH2*v25_KcF*v25_Kq1/v25_k8 + Cytc3p*v25_KmA*v25_Kq2 + QH2*v25_Kb1*v25_KcF*v25_Kq1/v25_k8 + v25_Kb2*v25_KmA*v25_Kq2) + Cytc3p*QH2 + Cytc3p*v25_KmA + QH2*v25_KmB);
            break;
        case 127:
            dwdp[18] = -1.0*Cytc2p*Cytc3p*QH2*v25_KcF*v25_v25_Complex_III*(Cytc3p*QH2*v25_Kq1/v25_k8 + QH2*v25_Kb1*v25_Kq1/v25_k8)/pow(Cytc2p*(Cytc3p*QH2*v25_KcF*v25_Kq1/v25_k8 + Cytc3p*v25_KmA*v25_Kq2 + QH2*v25_Kb1*v25_KcF*v25_Kq1/v25_k8 + v25_Kb2*v25_KmA*v25_Kq2) + Cytc3p*QH2 + Cytc3p*v25_KmA + QH2*v25_KmB, 2) + 1.0*Cytc3p*QH2*v25_v25_Complex_III/(Cytc2p*(Cytc3p*QH2*v25_KcF*v25_Kq1/v25_k8 + Cytc3p*v25_KmA*v25_Kq2 + QH2*v25_Kb1*v25_KcF*v25_Kq1/v25_k8 + v25_Kb2*v25_KmA*v25_Kq2) + Cytc3p*QH2 + Cytc3p*v25_KmA + QH2*v25_KmB);
            break;
        case 128:
            dwdp[18] = -1.0*Cytc2p*Cytc3p*QH2*v25_KcF*v25_v25_Complex_III*(-Cytc3p*QH2*v25_KcF*v25_Kq1/pow(v25_k8, 2) - QH2*v25_Kb1*v25_KcF*v25_Kq1/pow(v25_k8, 2))/pow(Cytc2p*(Cytc3p*QH2*v25_KcF*v25_Kq1/v25_k8 + Cytc3p*v25_KmA*v25_Kq2 + QH2*v25_Kb1*v25_KcF*v25_Kq1/v25_k8 + v25_Kb2*v25_KmA*v25_Kq2) + Cytc3p*QH2 + Cytc3p*v25_KmA + QH2*v25_KmB, 2);
            break;
        case 129:
            dwdp[18] = -1.0*Cytc2p*Cytc3p*QH2*v25_KcF*v25_v25_Complex_III*(Cytc3p*v25_KmA + v25_Kb2*v25_KmA)/pow(Cytc2p*(Cytc3p*QH2*v25_KcF*v25_Kq1/v25_k8 + Cytc3p*v25_KmA*v25_Kq2 + QH2*v25_Kb1*v25_KcF*v25_Kq1/v25_k8 + v25_Kb2*v25_KmA*v25_Kq2) + Cytc3p*QH2 + Cytc3p*v25_KmA + QH2*v25_KmB, 2);
            break;
        case 130:
            dwdp[18] = -1.0*Cytc2p*Cytc3p*QH2*v25_KcF*v25_v25_Complex_III*(Cytc3p*QH2*v25_KcF/v25_k8 + QH2*v25_Kb1*v25_KcF/v25_k8)/pow(Cytc2p*(Cytc3p*QH2*v25_KcF*v25_Kq1/v25_k8 + Cytc3p*v25_KmA*v25_Kq2 + QH2*v25_Kb1*v25_KcF*v25_Kq1/v25_k8 + v25_Kb2*v25_KmA*v25_Kq2) + Cytc3p*QH2 + Cytc3p*v25_KmA + QH2*v25_KmB, 2);
            break;
        case 131:
            dwdp[18] = -1.0*Cytc2p*Cytc3p*QH2*v25_KcF*v25_KmA*v25_Kq2*v25_v25_Complex_III/pow(Cytc2p*(Cytc3p*QH2*v25_KcF*v25_Kq1/v25_k8 + Cytc3p*v25_KmA*v25_Kq2 + QH2*v25_Kb1*v25_KcF*v25_Kq1/v25_k8 + v25_Kb2*v25_KmA*v25_Kq2) + Cytc3p*QH2 + Cytc3p*v25_KmA + QH2*v25_KmB, 2);
            break;
        case 132:
            dwdp[18] = -1.0*Cytc2p*Cytc3p*pow(QH2, 2)*pow(v25_KcF, 2)*v25_Kq1*v25_v25_Complex_III/(v25_k8*pow(Cytc2p*(Cytc3p*QH2*v25_KcF*v25_Kq1/v25_k8 + Cytc3p*v25_KmA*v25_Kq2 + QH2*v25_Kb1*v25_KcF*v25_Kq1/v25_k8 + v25_Kb2*v25_KmA*v25_Kq2) + Cytc3p*QH2 + Cytc3p*v25_KmA + QH2*v25_KmB, 2));
            break;
        case 133:
            dwdp[18] = -1.0*Cytc3p*pow(QH2, 2)*v25_KcF*v25_v25_Complex_III/pow(Cytc2p*(Cytc3p*QH2*v25_KcF*v25_Kq1/v25_k8 + Cytc3p*v25_KmA*v25_Kq2 + QH2*v25_Kb1*v25_KcF*v25_Kq1/v25_k8 + v25_Kb2*v25_KmA*v25_Kq2) + Cytc3p*QH2 + Cytc3p*v25_KmA + QH2*v25_KmB, 2);
            break;
        case 134:
            dwdp[18] = 1.0*Cytc3p*QH2*v25_KcF*v25_v25_Complex_III*(-Cytc2p*(Cytc3p*v25_Kq2 + v25_Kb2*v25_Kq2) - Cytc3p)/pow(Cytc2p*(Cytc3p*QH2*v25_KcF*v25_Kq1/v25_k8 + Cytc3p*v25_KmA*v25_Kq2 + QH2*v25_Kb1*v25_KcF*v25_Kq1/v25_k8 + v25_Kb2*v25_KmA*v25_Kq2) + Cytc3p*QH2 + Cytc3p*v25_KmA + QH2*v25_KmB, 2);
            break;
        case 135:
            dwdp[19] = 1.0*Cytc2p*v26_KcF/(Cytc2p + v26_Ks);
            break;
        case 136:
            dwdp[19] = 1.0*Cytc2p*v26_v26_Complex_IV/(Cytc2p + v26_Ks);
            break;
        case 137:
            dwdp[19] = -1.0*Cytc2p*v26_KcF*v26_v26_Complex_IV/pow(Cytc2p + v26_Ks, 2);
            break;
        case 138:
            dwdp[20] = 1.0*Acetyl_CoA_cyt*OXA_cyt*v27_Kc*v27_Kid*v27_V/(v27_Kb*v27_Keq*v27_Kia*(Acetyl_CoA_cyt*OXA_cyt + Acetyl_CoA_cyt*v27_Kb + OXA_cyt*v27_Ka + v27_Kia*v27_Kib));
            break;
        case 139:
            dwdp[20] = -1.0*Acetyl_CoA_cyt*OXA_cyt*v27_Kc*v27_Kid*v27_V*v27_v10_CS/(v27_Kb*pow(v27_Keq, 2)*v27_Kia*(Acetyl_CoA_cyt*OXA_cyt + Acetyl_CoA_cyt*v27_Kb + OXA_cyt*v27_Ka + v27_Kia*v27_Kib));
            break;
        case 140:
            dwdp[20] = 1.0*Acetyl_CoA_cyt*OXA_cyt*v27_Kc*v27_V*v27_v10_CS/(v27_Kb*v27_Keq*v27_Kia*(Acetyl_CoA_cyt*OXA_cyt + Acetyl_CoA_cyt*v27_Kb + OXA_cyt*v27_Ka + v27_Kia*v27_Kib));
            break;
        case 141:
            dwdp[20] = -1.0*Acetyl_CoA_cyt*OXA_cyt*v27_Kc*v27_Kid*v27_V*v27_v10_CS/(v27_Kb*v27_Keq*pow(Acetyl_CoA_cyt*OXA_cyt + Acetyl_CoA_cyt*v27_Kb + OXA_cyt*v27_Ka + v27_Kia*v27_Kib, 2));
            break;
        case 142:
            dwdp[20] = -1.0*Acetyl_CoA_cyt*OXA_cyt*v27_Kc*v27_Kib*v27_Kid*v27_V*v27_v10_CS/(v27_Kb*v27_Keq*v27_Kia*pow(Acetyl_CoA_cyt*OXA_cyt + Acetyl_CoA_cyt*v27_Kb + OXA_cyt*v27_Ka + v27_Kia*v27_Kib, 2)) - 1.0*Acetyl_CoA_cyt*OXA_cyt*v27_Kc*v27_Kid*v27_V*v27_v10_CS/(v27_Kb*v27_Keq*pow(v27_Kia, 2)*(Acetyl_CoA_cyt*OXA_cyt + Acetyl_CoA_cyt*v27_Kb + OXA_cyt*v27_Ka + v27_Kia*v27_Kib));
            break;
        case 143:
            dwdp[20] = 1.0*Acetyl_CoA_cyt*OXA_cyt*v27_Kid*v27_V*v27_v10_CS/(v27_Kb*v27_Keq*v27_Kia*(Acetyl_CoA_cyt*OXA_cyt + Acetyl_CoA_cyt*v27_Kb + OXA_cyt*v27_Ka + v27_Kia*v27_Kib));
            break;
        case 144:
            dwdp[20] = -1.0*pow(Acetyl_CoA_cyt, 2)*OXA_cyt*v27_Kc*v27_Kid*v27_V*v27_v10_CS/(v27_Kb*v27_Keq*v27_Kia*pow(Acetyl_CoA_cyt*OXA_cyt + Acetyl_CoA_cyt*v27_Kb + OXA_cyt*v27_Ka + v27_Kia*v27_Kib, 2)) - 1.0*Acetyl_CoA_cyt*OXA_cyt*v27_Kc*v27_Kid*v27_V*v27_v10_CS/(pow(v27_Kb, 2)*v27_Keq*v27_Kia*(Acetyl_CoA_cyt*OXA_cyt + Acetyl_CoA_cyt*v27_Kb + OXA_cyt*v27_Ka + v27_Kia*v27_Kib));
            break;
        case 145:
            dwdp[20] = -1.0*Acetyl_CoA_cyt*pow(OXA_cyt, 2)*v27_Kc*v27_Kid*v27_V*v27_v10_CS/(v27_Kb*v27_Keq*v27_Kia*pow(Acetyl_CoA_cyt*OXA_cyt + Acetyl_CoA_cyt*v27_Kb + OXA_cyt*v27_Ka + v27_Kia*v27_Kib, 2));
            break;
        case 146:
            dwdp[20] = 1.0*Acetyl_CoA_cyt*OXA_cyt*v27_Kc*v27_Kid*v27_v10_CS/(v27_Kb*v27_Keq*v27_Kia*(Acetyl_CoA_cyt*OXA_cyt + Acetyl_CoA_cyt*v27_Kb + OXA_cyt*v27_Ka + v27_Kia*v27_Kib));
            break;
        case 147:
            dwdp[21] = 1.0*ADP*v28_V/(pow(ADP, 2)/v28_Ki + ADP + v28_Km);
            break;
        case 148:
            dwdp[21] = 1.0*pow(ADP, 3)*v28_V*v28_v28_Complex_V/(pow(v28_Ki, 2)*pow(pow(ADP, 2)/v28_Ki + ADP + v28_Km, 2));
            break;
        case 149:
            dwdp[21] = -1.0*ADP*v28_V*v28_v28_Complex_V/pow(pow(ADP, 2)/v28_Ki + ADP + v28_Km, 2);
            break;
        case 150:
            dwdp[21] = 1.0*ADP*v28_v28_Complex_V/(pow(ADP, 2)/v28_Ki + ADP + v28_Km);
            break;
        case 151:
            dwdp[22] = 1.0*(Cit_cyt*v29_KcF*v29_Kp - IsoCitcyt*v29_KcR*v29_Ks)/(Cit_cyt*v29_Kp + IsoCitcyt*v29_Ks + v29_Kp*v29_Ks);
            break;
        case 152:
            dwdp[22] = -1.0*IsoCitcyt*v29_Ks*v29_v29_ACO/(Cit_cyt*v29_Kp + IsoCitcyt*v29_Ks + v29_Kp*v29_Ks);
            break;
        case 153:
            dwdp[22] = 1.0*Cit_cyt*v29_Kp*v29_v29_ACO/(Cit_cyt*v29_Kp + IsoCitcyt*v29_Ks + v29_Kp*v29_Ks);
            break;
        case 154:
            dwdp[22] = 1.0*Cit_cyt*v29_KcF*v29_v29_ACO/(Cit_cyt*v29_Kp + IsoCitcyt*v29_Ks + v29_Kp*v29_Ks) + 1.0*v29_v29_ACO*(-Cit_cyt - v29_Ks)*(Cit_cyt*v29_KcF*v29_Kp - IsoCitcyt*v29_KcR*v29_Ks)/pow(Cit_cyt*v29_Kp + IsoCitcyt*v29_Ks + v29_Kp*v29_Ks, 2);
            break;
        case 155:
            dwdp[22] = -1.0*IsoCitcyt*v29_KcR*v29_v29_ACO/(Cit_cyt*v29_Kp + IsoCitcyt*v29_Ks + v29_Kp*v29_Ks) + 1.0*v29_v29_ACO*(-IsoCitcyt - v29_Kp)*(Cit_cyt*v29_KcF*v29_Kp - IsoCitcyt*v29_KcR*v29_Ks)/pow(Cit_cyt*v29_Kp + IsoCitcyt*v29_Ks + v29_Kp*v29_Ks, 2);
            break;
        case 156:
            dwdp[23] = -1.0*pow(GAP, 2);
            break;
        case 157:
            dwdp[23] = 1.0*FBP;
            break;
        case 158:
            dwdp[24] = 1.0*(-Mal*OG_cyt*v30_KcR/(v30_KiP1*v30_KiP2*v30_beta) + Mal_cyt*OG*v30_KcF/(v30_KiS1*v30_KiS2*v30_alpha))/(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1);
            break;
        case 159:
            dwdp[24] = 1.0*Mal*OG*v30_v30_OGC*(-Mal*OG_cyt*v30_KcR/(v30_KiP1*v30_KiP2*v30_beta) + Mal_cyt*OG*v30_KcF/(v30_KiS1*v30_KiS2*v30_alpha))/(v30_KiP1*v30_KiS1*pow(v30_delta, 2)*pow(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1, 2));
            break;
        case 160:
            dwdp[24] = 1.0*Mal_cyt*OG_cyt*v30_v30_OGC*(-Mal*OG_cyt*v30_KcR/(v30_KiP1*v30_KiP2*v30_beta) + Mal_cyt*OG*v30_KcF/(v30_KiS1*v30_KiS2*v30_alpha))/(v30_KiP2*v30_KiS2*pow(v30_gamma, 2)*pow(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1, 2));
            break;
        case 161:
            dwdp[24] = 1.0*Mal*OG_cyt*v30_KcR*v30_v30_OGC/(v30_KiP1*v30_KiP2*pow(v30_beta, 2)*(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1)) + 1.0*Mal*OG_cyt*v30_v30_OGC*(-Mal*OG_cyt*v30_KcR/(v30_KiP1*v30_KiP2*v30_beta) + Mal_cyt*OG*v30_KcF/(v30_KiS1*v30_KiS2*v30_alpha))/(v30_KiP1*v30_KiP2*pow(v30_beta, 2)*pow(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1, 2));
            break;
        case 162:
            dwdp[24] = -1.0*Mal_cyt*OG*v30_KcF*v30_v30_OGC/(v30_KiS1*v30_KiS2*pow(v30_alpha, 2)*(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1)) + 1.0*Mal_cyt*OG*v30_v30_OGC*(-Mal*OG_cyt*v30_KcR/(v30_KiP1*v30_KiP2*v30_beta) + Mal_cyt*OG*v30_KcF/(v30_KiS1*v30_KiS2*v30_alpha))/(v30_KiS1*v30_KiS2*pow(v30_alpha, 2)*pow(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1, 2));
            break;
        case 163:
            dwdp[24] = -1.0*Mal*OG_cyt*v30_v30_OGC/(v30_KiP1*v30_KiP2*v30_beta*(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1));
            break;
        case 164:
            dwdp[24] = 1.0*Mal_cyt*OG*v30_v30_OGC/(v30_KiS1*v30_KiS2*v30_alpha*(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1));
            break;
        case 165:
            dwdp[24] = 1.0*Mal*OG_cyt*v30_KcR*v30_v30_OGC/(v30_KiP1*pow(v30_KiP2, 2)*v30_beta*(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1)) + 1.0*v30_v30_OGC*(-Mal*OG_cyt*v30_KcR/(v30_KiP1*v30_KiP2*v30_beta) + Mal_cyt*OG*v30_KcF/(v30_KiS1*v30_KiS2*v30_alpha))*(Mal*OG_cyt/(v30_KiP1*pow(v30_KiP2, 2)*v30_beta) + Mal_cyt*OG_cyt/(pow(v30_KiP2, 2)*v30_KiS2*v30_gamma) + OG_cyt/pow(v30_KiP2, 2))/pow(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1, 2);
            break;
        case 166:
            dwdp[24] = 1.0*Mal*OG_cyt*v30_KcR*v30_v30_OGC/(pow(v30_KiP1, 2)*v30_KiP2*v30_beta*(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1)) + 1.0*v30_v30_OGC*(-Mal*OG_cyt*v30_KcR/(v30_KiP1*v30_KiP2*v30_beta) + Mal_cyt*OG*v30_KcF/(v30_KiS1*v30_KiS2*v30_alpha))*(Mal*OG/(pow(v30_KiP1, 2)*v30_KiS1*v30_delta) + Mal*OG_cyt/(pow(v30_KiP1, 2)*v30_KiP2*v30_beta) + Mal/pow(v30_KiP1, 2))/pow(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1, 2);
            break;
        case 167:
            dwdp[24] = -1.0*Mal_cyt*OG*v30_KcF*v30_v30_OGC/(v30_KiS1*pow(v30_KiS2, 2)*v30_alpha*(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1)) + 1.0*v30_v30_OGC*(-Mal*OG_cyt*v30_KcR/(v30_KiP1*v30_KiP2*v30_beta) + Mal_cyt*OG*v30_KcF/(v30_KiS1*v30_KiS2*v30_alpha))*(Mal_cyt*OG/(v30_KiS1*pow(v30_KiS2, 2)*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*pow(v30_KiS2, 2)*v30_gamma) + Mal_cyt/pow(v30_KiS2, 2))/pow(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1, 2);
            break;
        case 168:
            dwdp[24] = -1.0*Mal_cyt*OG*v30_KcF*v30_v30_OGC/(pow(v30_KiS1, 2)*v30_KiS2*v30_alpha*(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1)) + 1.0*v30_v30_OGC*(-Mal*OG_cyt*v30_KcR/(v30_KiP1*v30_KiP2*v30_beta) + Mal_cyt*OG*v30_KcF/(v30_KiS1*v30_KiS2*v30_alpha))*(Mal*OG/(v30_KiP1*pow(v30_KiS1, 2)*v30_delta) + Mal_cyt*OG/(pow(v30_KiS1, 2)*v30_KiS2*v30_alpha) + OG/pow(v30_KiS1, 2))/pow(Mal*OG/(v30_KiP1*v30_KiS1*v30_delta) + Mal*OG_cyt/(v30_KiP1*v30_KiP2*v30_beta) + Mal/v30_KiP1 + Mal_cyt*OG/(v30_KiS1*v30_KiS2*v30_alpha) + Mal_cyt*OG_cyt/(v30_KiP2*v30_KiS2*v30_gamma) + Mal_cyt/v30_KiS2 + OG/v30_KiS1 + OG_cyt/v30_KiP2 + 1, 2);
            break;
        case 169:
            dwdp[25] = 1.0*(-Mal_cyt*NAD*v31_kminus1*v31_kminus2*v31_kminus3*v31_kminus4 + NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_k3*v31_k4)/(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2);
            break;
        case 170:
            dwdp[25] = -1.0*Mal_cyt*NAD*v31_kminus1*v31_kminus2*v31_kminus3*v31_v31_MDH/(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2) + 1.0*v31_v31_MDH*(-Mal_cyt*NAD*v31_kminus1*v31_kminus2*v31_kminus3*v31_kminus4 + NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_k3*v31_k4)*(-Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3 - Mal_cyt*NAD*v31_kminus3*(v31_kminus1 + v31_kminus2) - NAD*OXA_cyt*v31_k2*v31_k3 - NAD*v31_kminus1*(v31_k3 + v31_kminus2))/pow(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2, 2);
            break;
        case 171:
            dwdp[25] = -1.0*Mal_cyt*NAD*v31_kminus1*v31_kminus2*v31_kminus4*v31_v31_MDH/(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2) + 1.0*v31_v31_MDH*(-Mal_cyt*NAD*v31_kminus1*v31_kminus2*v31_kminus3*v31_kminus4 + NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_k3*v31_k4)*(-Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus4 - Mal_cyt*NAD*v31_kminus4*(v31_kminus1 + v31_kminus2) - Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2 - Mal_cyt*NADH_cyt - Mal_cyt*v31_kminus1*v31_kminus2)/pow(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2, 2);
            break;
        case 172:
            dwdp[25] = -1.0*Mal_cyt*NAD*v31_kminus1*v31_kminus3*v31_kminus4*v31_v31_MDH/(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2) + 1.0*v31_v31_MDH*(-Mal_cyt*NAD*v31_kminus1*v31_kminus2*v31_kminus3*v31_kminus4 + NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_k3*v31_k4)*(-Mal_cyt*NAD*v31_kminus3*v31_kminus4 - Mal_cyt*v31_kminus1*v31_kminus3 - NAD*v31_kminus1*v31_kminus4 - NADH_cyt*v31_k1*v31_k4 - v31_k4*v31_kminus1 - 1)/pow(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2, 2);
            break;
        case 173:
            dwdp[25] = -1.0*Mal_cyt*NAD*v31_kminus2*v31_kminus3*v31_kminus4*v31_v31_MDH/(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2) + 1.0*v31_v31_MDH*(-Mal_cyt*NAD*v31_kminus1*v31_kminus2*v31_kminus3*v31_kminus4 + NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_k3*v31_k4)*(-Mal_cyt*NAD*v31_kminus3*v31_kminus4 - Mal_cyt*v31_kminus2*v31_kminus3 - NAD*v31_kminus4*(v31_k3 + v31_kminus2) - v31_k4*(v31_k3 + v31_kminus2))/pow(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2, 2);
            break;
        case 174:
            dwdp[25] = 1.0*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_k3*v31_v31_MDH/(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2) + 1.0*v31_v31_MDH*(-Mal_cyt*NAD*v31_kminus1*v31_kminus2*v31_kminus3*v31_kminus4 + NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_k3*v31_k4)*(-NADH_cyt*OXA_cyt*v31_k1*v31_k2 - NADH_cyt*v31_k1*(v31_k3 + v31_kminus2) - OXA_cyt*v31_k2*v31_k3 - v31_kminus1*(v31_k3 + v31_kminus2))/pow(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2, 2);
            break;
        case 175:
            dwdp[25] = 1.0*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_k4*v31_v31_MDH/(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2) + 1.0*v31_v31_MDH*(-Mal_cyt*NAD*v31_kminus1*v31_kminus2*v31_kminus3*v31_kminus4 + NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_k3*v31_k4)*(-NAD*OXA_cyt*v31_k2*v31_kminus4 - NAD*v31_kminus1*v31_kminus4 - NADH_cyt*OXA_cyt*v31_k1*v31_k2 - NADH_cyt*v31_k1*v31_k4 - OXA_cyt*v31_k2*v31_k4 - v31_k4*v31_kminus1)/pow(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2, 2);
            break;
        case 176:
            dwdp[25] = 1.0*NADH_cyt*OXA_cyt*v31_k1*v31_k3*v31_k4*v31_v31_MDH/(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2) + 1.0*v31_v31_MDH*(-Mal_cyt*NAD*v31_kminus1*v31_kminus2*v31_kminus3*v31_kminus4 + NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_k3*v31_k4)*(-Mal_cyt*NAD*OXA_cyt*v31_kminus3*v31_kminus4 - Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_kminus3 - NAD*OXA_cyt*v31_k3*v31_kminus4 - NADH_cyt*OXA_cyt*v31_k1*(v31_k3 + v31_k4) - OXA_cyt*v31_k3*v31_k4)/pow(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2, 2);
            break;
        case 177:
            dwdp[25] = 1.0*NADH_cyt*OXA_cyt*v31_k2*v31_k3*v31_k4*v31_v31_MDH/(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2) + 1.0*v31_v31_MDH*(-Mal_cyt*NAD*v31_kminus1*v31_kminus2*v31_kminus3*v31_kminus4 + NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_k3*v31_k4)*(-Mal_cyt*NADH_cyt*OXA_cyt*v31_k2*v31_kminus3 - NADH_cyt*OXA_cyt*v31_k2*(v31_k3 + v31_k4) - NADH_cyt*v31_k4*(v31_k3 + v31_kminus2) - 1)/pow(Mal_cyt*NAD*OXA_cyt*v31_k2*v31_kminus3*v31_kminus4 + Mal_cyt*NAD*v31_kminus3*v31_kminus4*(v31_kminus1 + v31_kminus2) + Mal_cyt*NADH_cyt*OXA_cyt*v31_k1*v31_k2*v31_kminus3 + Mal_cyt*NADH_cyt*v31_kminus3 + Mal_cyt*v31_kminus1*v31_kminus2*v31_kminus3 + NAD*OXA_cyt*v31_k2*v31_k3*v31_kminus4 + NAD*v31_kminus1*v31_kminus4*(v31_k3 + v31_kminus2) + NADH_cyt*OXA_cyt*v31_k1*v31_k2*(v31_k3 + v31_k4) + NADH_cyt*v31_k1*v31_k4*(v31_k3 + v31_kminus2) + OXA_cyt*v31_k2*v31_k3*v31_k4 + v31_k1 + v31_k4*v31_kminus1*(v31_k3 + v31_kminus2) + v31_kminus2, 2);
            break;
        case 178:
            dwdp[26] = 1.0*v32_KcF*v32_KcR*(Asp_cyt*OG_cyt - Glu_cyt*OXA_cyt/v32_Keq)/(Asp_cyt*OG_cyt*v32_KcR + Asp_cyt*OXA_cyt*v32_KcF*v32_KmP2/(v32_Keq*v32_KiS1) + Asp_cyt*v32_KcR*v32_KmS2 + Glu_cyt*OG_cyt*v32_KcR*v32_KmS1/v32_KiP2 + Glu_cyt*OXA_cyt*v32_KcF/v32_Keq + Glu_cyt*v32_KcF*v32_KmP1/v32_Keq + OG_cyt*v32_KcR*v32_KmS1 + OXA_cyt*v32_KcF*v32_KmP2/v32_Keq);
            break;
        case 179:
            dwdp[26] = 1.0*v32_KcF*v32_KcR*v32_v32_AspTA*(Asp_cyt*OG_cyt - Glu_cyt*OXA_cyt/v32_Keq)*(-Asp_cyt*OG_cyt - Asp_cyt*v32_KmS2 - Glu_cyt*OG_cyt*v32_KmS1/v32_KiP2 - OG_cyt*v32_KmS1)/pow(Asp_cyt*OG_cyt*v32_KcR + Asp_cyt*OXA_cyt*v32_KcF*v32_KmP2/(v32_Keq*v32_KiS1) + Asp_cyt*v32_KcR*v32_KmS2 + Glu_cyt*OG_cyt*v32_KcR*v32_KmS1/v32_KiP2 + Glu_cyt*OXA_cyt*v32_KcF/v32_Keq + Glu_cyt*v32_KcF*v32_KmP1/v32_Keq + OG_cyt*v32_KcR*v32_KmS1 + OXA_cyt*v32_KcF*v32_KmP2/v32_Keq, 2) + 1.0*v32_KcF*v32_v32_AspTA*(Asp_cyt*OG_cyt - Glu_cyt*OXA_cyt/v32_Keq)/(Asp_cyt*OG_cyt*v32_KcR + Asp_cyt*OXA_cyt*v32_KcF*v32_KmP2/(v32_Keq*v32_KiS1) + Asp_cyt*v32_KcR*v32_KmS2 + Glu_cyt*OG_cyt*v32_KcR*v32_KmS1/v32_KiP2 + Glu_cyt*OXA_cyt*v32_KcF/v32_Keq + Glu_cyt*v32_KcF*v32_KmP1/v32_Keq + OG_cyt*v32_KcR*v32_KmS1 + OXA_cyt*v32_KcF*v32_KmP2/v32_Keq);
            break;
        case 180:
            dwdp[26] = 1.0*v32_KcF*v32_KcR*v32_v32_AspTA*(Asp_cyt*OG_cyt - Glu_cyt*OXA_cyt/v32_Keq)*(-Asp_cyt*OXA_cyt*v32_KmP2/(v32_Keq*v32_KiS1) - Glu_cyt*OXA_cyt/v32_Keq - Glu_cyt*v32_KmP1/v32_Keq - OXA_cyt*v32_KmP2/v32_Keq)/pow(Asp_cyt*OG_cyt*v32_KcR + Asp_cyt*OXA_cyt*v32_KcF*v32_KmP2/(v32_Keq*v32_KiS1) + Asp_cyt*v32_KcR*v32_KmS2 + Glu_cyt*OG_cyt*v32_KcR*v32_KmS1/v32_KiP2 + Glu_cyt*OXA_cyt*v32_KcF/v32_Keq + Glu_cyt*v32_KcF*v32_KmP1/v32_Keq + OG_cyt*v32_KcR*v32_KmS1 + OXA_cyt*v32_KcF*v32_KmP2/v32_Keq, 2) + 1.0*v32_KcR*v32_v32_AspTA*(Asp_cyt*OG_cyt - Glu_cyt*OXA_cyt/v32_Keq)/(Asp_cyt*OG_cyt*v32_KcR + Asp_cyt*OXA_cyt*v32_KcF*v32_KmP2/(v32_Keq*v32_KiS1) + Asp_cyt*v32_KcR*v32_KmS2 + Glu_cyt*OG_cyt*v32_KcR*v32_KmS1/v32_KiP2 + Glu_cyt*OXA_cyt*v32_KcF/v32_Keq + Glu_cyt*v32_KcF*v32_KmP1/v32_Keq + OG_cyt*v32_KcR*v32_KmS1 + OXA_cyt*v32_KcF*v32_KmP2/v32_Keq);
            break;
        case 181:
            dwdp[26] = 1.0*Glu_cyt*OXA_cyt*v32_KcF*v32_KcR*v32_v32_AspTA/(pow(v32_Keq, 2)*(Asp_cyt*OG_cyt*v32_KcR + Asp_cyt*OXA_cyt*v32_KcF*v32_KmP2/(v32_Keq*v32_KiS1) + Asp_cyt*v32_KcR*v32_KmS2 + Glu_cyt*OG_cyt*v32_KcR*v32_KmS1/v32_KiP2 + Glu_cyt*OXA_cyt*v32_KcF/v32_Keq + Glu_cyt*v32_KcF*v32_KmP1/v32_Keq + OG_cyt*v32_KcR*v32_KmS1 + OXA_cyt*v32_KcF*v32_KmP2/v32_Keq)) + 1.0*v32_KcF*v32_KcR*v32_v32_AspTA*(Asp_cyt*OG_cyt - Glu_cyt*OXA_cyt/v32_Keq)*(Asp_cyt*OXA_cyt*v32_KcF*v32_KmP2/(pow(v32_Keq, 2)*v32_KiS1) + Glu_cyt*OXA_cyt*v32_KcF/pow(v32_Keq, 2) + Glu_cyt*v32_KcF*v32_KmP1/pow(v32_Keq, 2) + OXA_cyt*v32_KcF*v32_KmP2/pow(v32_Keq, 2))/pow(Asp_cyt*OG_cyt*v32_KcR + Asp_cyt*OXA_cyt*v32_KcF*v32_KmP2/(v32_Keq*v32_KiS1) + Asp_cyt*v32_KcR*v32_KmS2 + Glu_cyt*OG_cyt*v32_KcR*v32_KmS1/v32_KiP2 + Glu_cyt*OXA_cyt*v32_KcF/v32_Keq + Glu_cyt*v32_KcF*v32_KmP1/v32_Keq + OG_cyt*v32_KcR*v32_KmS1 + OXA_cyt*v32_KcF*v32_KmP2/v32_Keq, 2);
            break;
        case 182:
            dwdp[26] = 1.0*Glu_cyt*OG_cyt*v32_KcF*pow(v32_KcR, 2)*v32_KmS1*v32_v32_AspTA*(Asp_cyt*OG_cyt - Glu_cyt*OXA_cyt/v32_Keq)/(pow(v32_KiP2, 2)*pow(Asp_cyt*OG_cyt*v32_KcR + Asp_cyt*OXA_cyt*v32_KcF*v32_KmP2/(v32_Keq*v32_KiS1) + Asp_cyt*v32_KcR*v32_KmS2 + Glu_cyt*OG_cyt*v32_KcR*v32_KmS1/v32_KiP2 + Glu_cyt*OXA_cyt*v32_KcF/v32_Keq + Glu_cyt*v32_KcF*v32_KmP1/v32_Keq + OG_cyt*v32_KcR*v32_KmS1 + OXA_cyt*v32_KcF*v32_KmP2/v32_Keq, 2));
            break;
        case 183:
            dwdp[26] = 1.0*Asp_cyt*OXA_cyt*pow(v32_KcF, 2)*v32_KcR*v32_KmP2*v32_v32_AspTA*(Asp_cyt*OG_cyt - Glu_cyt*OXA_cyt/v32_Keq)/(v32_Keq*pow(v32_KiS1, 2)*pow(Asp_cyt*OG_cyt*v32_KcR + Asp_cyt*OXA_cyt*v32_KcF*v32_KmP2/(v32_Keq*v32_KiS1) + Asp_cyt*v32_KcR*v32_KmS2 + Glu_cyt*OG_cyt*v32_KcR*v32_KmS1/v32_KiP2 + Glu_cyt*OXA_cyt*v32_KcF/v32_Keq + Glu_cyt*v32_KcF*v32_KmP1/v32_Keq + OG_cyt*v32_KcR*v32_KmS1 + OXA_cyt*v32_KcF*v32_KmP2/v32_Keq, 2));
            break;
        case 184:
            dwdp[26] = 1.0*v32_KcF*v32_KcR*v32_v32_AspTA*(Asp_cyt*OG_cyt - Glu_cyt*OXA_cyt/v32_Keq)*(-Asp_cyt*OXA_cyt*v32_KcF/(v32_Keq*v32_KiS1) - OXA_cyt*v32_KcF/v32_Keq)/pow(Asp_cyt*OG_cyt*v32_KcR + Asp_cyt*OXA_cyt*v32_KcF*v32_KmP2/(v32_Keq*v32_KiS1) + Asp_cyt*v32_KcR*v32_KmS2 + Glu_cyt*OG_cyt*v32_KcR*v32_KmS1/v32_KiP2 + Glu_cyt*OXA_cyt*v32_KcF/v32_Keq + Glu_cyt*v32_KcF*v32_KmP1/v32_Keq + OG_cyt*v32_KcR*v32_KmS1 + OXA_cyt*v32_KcF*v32_KmP2/v32_Keq, 2);
            break;
        case 185:
            dwdp[26] = -1.0*Glu_cyt*pow(v32_KcF, 2)*v32_KcR*v32_v32_AspTA*(Asp_cyt*OG_cyt - Glu_cyt*OXA_cyt/v32_Keq)/(v32_Keq*pow(Asp_cyt*OG_cyt*v32_KcR + Asp_cyt*OXA_cyt*v32_KcF*v32_KmP2/(v32_Keq*v32_KiS1) + Asp_cyt*v32_KcR*v32_KmS2 + Glu_cyt*OG_cyt*v32_KcR*v32_KmS1/v32_KiP2 + Glu_cyt*OXA_cyt*v32_KcF/v32_Keq + Glu_cyt*v32_KcF*v32_KmP1/v32_Keq + OG_cyt*v32_KcR*v32_KmS1 + OXA_cyt*v32_KcF*v32_KmP2/v32_Keq, 2));
            break;
        case 186:
            dwdp[26] = -1.0*Asp_cyt*v32_KcF*pow(v32_KcR, 2)*v32_v32_AspTA*(Asp_cyt*OG_cyt - Glu_cyt*OXA_cyt/v32_Keq)/pow(Asp_cyt*OG_cyt*v32_KcR + Asp_cyt*OXA_cyt*v32_KcF*v32_KmP2/(v32_Keq*v32_KiS1) + Asp_cyt*v32_KcR*v32_KmS2 + Glu_cyt*OG_cyt*v32_KcR*v32_KmS1/v32_KiP2 + Glu_cyt*OXA_cyt*v32_KcF/v32_Keq + Glu_cyt*v32_KcF*v32_KmP1/v32_Keq + OG_cyt*v32_KcR*v32_KmS1 + OXA_cyt*v32_KcF*v32_KmP2/v32_Keq, 2);
            break;
        case 187:
            dwdp[26] = 1.0*v32_KcF*v32_KcR*v32_v32_AspTA*(Asp_cyt*OG_cyt - Glu_cyt*OXA_cyt/v32_Keq)*(-Glu_cyt*OG_cyt*v32_KcR/v32_KiP2 - OG_cyt*v32_KcR)/pow(Asp_cyt*OG_cyt*v32_KcR + Asp_cyt*OXA_cyt*v32_KcF*v32_KmP2/(v32_Keq*v32_KiS1) + Asp_cyt*v32_KcR*v32_KmS2 + Glu_cyt*OG_cyt*v32_KcR*v32_KmS1/v32_KiP2 + Glu_cyt*OXA_cyt*v32_KcF/v32_Keq + Glu_cyt*v32_KcF*v32_KmP1/v32_Keq + OG_cyt*v32_KcR*v32_KmS1 + OXA_cyt*v32_KcF*v32_KmP2/v32_Keq, 2);
            break;
        case 188:
            dwdp[27] = 1.0*(-Cit*Mal_cyt*v33_KcR/(v33_KiP1*v33_KiP2*v33_beta) + Cit_cyt*Mal*v33_KcF/(v33_KiS1*v33_KiS2*v33_alpha))/(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1);
            break;
        case 189:
            dwdp[27] = 1.0*Cit_cyt*Mal_cyt*v33_v33_CIC*(-Cit*Mal_cyt*v33_KcR/(v33_KiP1*v33_KiP2*v33_beta) + Cit_cyt*Mal*v33_KcF/(v33_KiS1*v33_KiS2*v33_alpha))/(v33_KiP1*v33_KiS1*pow(v33_delta, 2)*pow(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1, 2));
            break;
        case 190:
            dwdp[27] = 1.0*Cit*Mal*v33_v33_CIC*(-Cit*Mal_cyt*v33_KcR/(v33_KiP1*v33_KiP2*v33_beta) + Cit_cyt*Mal*v33_KcF/(v33_KiS1*v33_KiS2*v33_alpha))/(v33_KiP2*v33_KiS2*pow(v33_gamma, 2)*pow(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1, 2));
            break;
        case 191:
            dwdp[27] = 1.0*Cit*Mal_cyt*v33_KcR*v33_v33_CIC/(v33_KiP1*v33_KiP2*pow(v33_beta, 2)*(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1)) + 1.0*Cit*Mal_cyt*v33_v33_CIC*(-Cit*Mal_cyt*v33_KcR/(v33_KiP1*v33_KiP2*v33_beta) + Cit_cyt*Mal*v33_KcF/(v33_KiS1*v33_KiS2*v33_alpha))/(v33_KiP1*v33_KiP2*pow(v33_beta, 2)*pow(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1, 2));
            break;
        case 192:
            dwdp[27] = -1.0*Cit_cyt*Mal*v33_KcF*v33_v33_CIC/(v33_KiS1*v33_KiS2*pow(v33_alpha, 2)*(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1)) + 1.0*Cit_cyt*Mal*v33_v33_CIC*(-Cit*Mal_cyt*v33_KcR/(v33_KiP1*v33_KiP2*v33_beta) + Cit_cyt*Mal*v33_KcF/(v33_KiS1*v33_KiS2*v33_alpha))/(v33_KiS1*v33_KiS2*pow(v33_alpha, 2)*pow(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1, 2));
            break;
        case 193:
            dwdp[27] = -1.0*Cit*Mal_cyt*v33_v33_CIC/(v33_KiP1*v33_KiP2*v33_beta*(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1));
            break;
        case 194:
            dwdp[27] = 1.0*Cit_cyt*Mal*v33_v33_CIC/(v33_KiS1*v33_KiS2*v33_alpha*(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1));
            break;
        case 195:
            dwdp[27] = 1.0*Cit*Mal_cyt*v33_KcR*v33_v33_CIC/(v33_KiP1*pow(v33_KiP2, 2)*v33_beta*(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1)) + 1.0*v33_v33_CIC*(-Cit*Mal_cyt*v33_KcR/(v33_KiP1*v33_KiP2*v33_beta) + Cit_cyt*Mal*v33_KcF/(v33_KiS1*v33_KiS2*v33_alpha))*(Cit*Mal/(pow(v33_KiP2, 2)*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*pow(v33_KiP2, 2)*v33_beta) + Cit/pow(v33_KiP2, 2))/pow(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1, 2);
            break;
        case 196:
            dwdp[27] = 1.0*Cit*Mal_cyt*v33_KcR*v33_v33_CIC/(pow(v33_KiP1, 2)*v33_KiP2*v33_beta*(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1)) + 1.0*v33_v33_CIC*(-Cit*Mal_cyt*v33_KcR/(v33_KiP1*v33_KiP2*v33_beta) + Cit_cyt*Mal*v33_KcF/(v33_KiS1*v33_KiS2*v33_alpha))*(Cit*Mal_cyt/(pow(v33_KiP1, 2)*v33_KiP2*v33_beta) + Cit_cyt*Mal_cyt/(pow(v33_KiP1, 2)*v33_KiS1*v33_delta) + Mal_cyt/pow(v33_KiP1, 2))/pow(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1, 2);
            break;
        case 197:
            dwdp[27] = -1.0*Cit_cyt*Mal*v33_KcF*v33_v33_CIC/(v33_KiS1*pow(v33_KiS2, 2)*v33_alpha*(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1)) + 1.0*v33_v33_CIC*(-Cit*Mal_cyt*v33_KcR/(v33_KiP1*v33_KiP2*v33_beta) + Cit_cyt*Mal*v33_KcF/(v33_KiS1*v33_KiS2*v33_alpha))*(Cit*Mal/(v33_KiP2*pow(v33_KiS2, 2)*v33_gamma) + Cit_cyt*Mal/(v33_KiS1*pow(v33_KiS2, 2)*v33_alpha) + Mal/pow(v33_KiS2, 2))/pow(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1, 2);
            break;
        case 198:
            dwdp[27] = -1.0*Cit_cyt*Mal*v33_KcF*v33_v33_CIC/(pow(v33_KiS1, 2)*v33_KiS2*v33_alpha*(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1)) + 1.0*v33_v33_CIC*(-Cit*Mal_cyt*v33_KcR/(v33_KiP1*v33_KiP2*v33_beta) + Cit_cyt*Mal*v33_KcF/(v33_KiS1*v33_KiS2*v33_alpha))*(Cit_cyt*Mal/(pow(v33_KiS1, 2)*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*pow(v33_KiS1, 2)*v33_delta) + Cit_cyt/pow(v33_KiS1, 2))/pow(Cit*Mal/(v33_KiP2*v33_KiS2*v33_gamma) + Cit*Mal_cyt/(v33_KiP1*v33_KiP2*v33_beta) + Cit/v33_KiP2 + Cit_cyt*Mal/(v33_KiS1*v33_KiS2*v33_alpha) + Cit_cyt*Mal_cyt/(v33_KiP1*v33_KiS1*v33_delta) + Cit_cyt/v33_KiS1 + Mal/v33_KiS2 + Mal_cyt/v33_KiP1 + 1, 2);
            break;
        case 199:
            dwdp[28] = 1.0*v34_KcF*v34_KcR*(-ETFox*QH2/v34_Keq + ETFred*Q)/(ETFox*ETFred*v34_KcF*v34_KmP2/(v34_Keq*v34_KiS1) + ETFox*QH2*v34_KcF/v34_Keq + ETFox*v34_KcF*v34_KmP2/v34_Keq + ETFred*Q*v34_KcR + ETFred*v34_KcR*v34_KmS2 + Q*QH2*v34_KcR*v34_KmS1/v34_KiP2 + Q*v34_KcR*v34_KmS1 + QH2*v34_KcF*v34_KmP1/v34_Keq);
            break;
        case 200:
            dwdp[28] = 1.0*v34_KcF*v34_KcR*v34_v34_ETF_QO*(-ETFox*QH2/v34_Keq + ETFred*Q)*(-ETFred*Q - ETFred*v34_KmS2 - Q*QH2*v34_KmS1/v34_KiP2 - Q*v34_KmS1)/pow(ETFox*ETFred*v34_KcF*v34_KmP2/(v34_Keq*v34_KiS1) + ETFox*QH2*v34_KcF/v34_Keq + ETFox*v34_KcF*v34_KmP2/v34_Keq + ETFred*Q*v34_KcR + ETFred*v34_KcR*v34_KmS2 + Q*QH2*v34_KcR*v34_KmS1/v34_KiP2 + Q*v34_KcR*v34_KmS1 + QH2*v34_KcF*v34_KmP1/v34_Keq, 2) + 1.0*v34_KcF*v34_v34_ETF_QO*(-ETFox*QH2/v34_Keq + ETFred*Q)/(ETFox*ETFred*v34_KcF*v34_KmP2/(v34_Keq*v34_KiS1) + ETFox*QH2*v34_KcF/v34_Keq + ETFox*v34_KcF*v34_KmP2/v34_Keq + ETFred*Q*v34_KcR + ETFred*v34_KcR*v34_KmS2 + Q*QH2*v34_KcR*v34_KmS1/v34_KiP2 + Q*v34_KcR*v34_KmS1 + QH2*v34_KcF*v34_KmP1/v34_Keq);
            break;
        case 201:
            dwdp[28] = 1.0*v34_KcF*v34_KcR*v34_v34_ETF_QO*(-ETFox*QH2/v34_Keq + ETFred*Q)*(-ETFox*ETFred*v34_KmP2/(v34_Keq*v34_KiS1) - ETFox*QH2/v34_Keq - ETFox*v34_KmP2/v34_Keq - QH2*v34_KmP1/v34_Keq)/pow(ETFox*ETFred*v34_KcF*v34_KmP2/(v34_Keq*v34_KiS1) + ETFox*QH2*v34_KcF/v34_Keq + ETFox*v34_KcF*v34_KmP2/v34_Keq + ETFred*Q*v34_KcR + ETFred*v34_KcR*v34_KmS2 + Q*QH2*v34_KcR*v34_KmS1/v34_KiP2 + Q*v34_KcR*v34_KmS1 + QH2*v34_KcF*v34_KmP1/v34_Keq, 2) + 1.0*v34_KcR*v34_v34_ETF_QO*(-ETFox*QH2/v34_Keq + ETFred*Q)/(ETFox*ETFred*v34_KcF*v34_KmP2/(v34_Keq*v34_KiS1) + ETFox*QH2*v34_KcF/v34_Keq + ETFox*v34_KcF*v34_KmP2/v34_Keq + ETFred*Q*v34_KcR + ETFred*v34_KcR*v34_KmS2 + Q*QH2*v34_KcR*v34_KmS1/v34_KiP2 + Q*v34_KcR*v34_KmS1 + QH2*v34_KcF*v34_KmP1/v34_Keq);
            break;
        case 202:
            dwdp[28] = 1.0*ETFox*QH2*v34_KcF*v34_KcR*v34_v34_ETF_QO/(pow(v34_Keq, 2)*(ETFox*ETFred*v34_KcF*v34_KmP2/(v34_Keq*v34_KiS1) + ETFox*QH2*v34_KcF/v34_Keq + ETFox*v34_KcF*v34_KmP2/v34_Keq + ETFred*Q*v34_KcR + ETFred*v34_KcR*v34_KmS2 + Q*QH2*v34_KcR*v34_KmS1/v34_KiP2 + Q*v34_KcR*v34_KmS1 + QH2*v34_KcF*v34_KmP1/v34_Keq)) + 1.0*v34_KcF*v34_KcR*v34_v34_ETF_QO*(-ETFox*QH2/v34_Keq + ETFred*Q)*(ETFox*ETFred*v34_KcF*v34_KmP2/(pow(v34_Keq, 2)*v34_KiS1) + ETFox*QH2*v34_KcF/pow(v34_Keq, 2) + ETFox*v34_KcF*v34_KmP2/pow(v34_Keq, 2) + QH2*v34_KcF*v34_KmP1/pow(v34_Keq, 2))/pow(ETFox*ETFred*v34_KcF*v34_KmP2/(v34_Keq*v34_KiS1) + ETFox*QH2*v34_KcF/v34_Keq + ETFox*v34_KcF*v34_KmP2/v34_Keq + ETFred*Q*v34_KcR + ETFred*v34_KcR*v34_KmS2 + Q*QH2*v34_KcR*v34_KmS1/v34_KiP2 + Q*v34_KcR*v34_KmS1 + QH2*v34_KcF*v34_KmP1/v34_Keq, 2);
            break;
        case 203:
            dwdp[28] = 1.0*Q*QH2*v34_KcF*pow(v34_KcR, 2)*v34_KmS1*v34_v34_ETF_QO*(-ETFox*QH2/v34_Keq + ETFred*Q)/(pow(v34_KiP2, 2)*pow(ETFox*ETFred*v34_KcF*v34_KmP2/(v34_Keq*v34_KiS1) + ETFox*QH2*v34_KcF/v34_Keq + ETFox*v34_KcF*v34_KmP2/v34_Keq + ETFred*Q*v34_KcR + ETFred*v34_KcR*v34_KmS2 + Q*QH2*v34_KcR*v34_KmS1/v34_KiP2 + Q*v34_KcR*v34_KmS1 + QH2*v34_KcF*v34_KmP1/v34_Keq, 2));
            break;
        case 204:
            dwdp[28] = 1.0*ETFox*ETFred*pow(v34_KcF, 2)*v34_KcR*v34_KmP2*v34_v34_ETF_QO*(-ETFox*QH2/v34_Keq + ETFred*Q)/(v34_Keq*pow(v34_KiS1, 2)*pow(ETFox*ETFred*v34_KcF*v34_KmP2/(v34_Keq*v34_KiS1) + ETFox*QH2*v34_KcF/v34_Keq + ETFox*v34_KcF*v34_KmP2/v34_Keq + ETFred*Q*v34_KcR + ETFred*v34_KcR*v34_KmS2 + Q*QH2*v34_KcR*v34_KmS1/v34_KiP2 + Q*v34_KcR*v34_KmS1 + QH2*v34_KcF*v34_KmP1/v34_Keq, 2));
            break;
        case 205:
            dwdp[28] = 1.0*v34_KcF*v34_KcR*v34_v34_ETF_QO*(-ETFox*QH2/v34_Keq + ETFred*Q)*(-ETFox*ETFred*v34_KcF/(v34_Keq*v34_KiS1) - ETFox*v34_KcF/v34_Keq)/pow(ETFox*ETFred*v34_KcF*v34_KmP2/(v34_Keq*v34_KiS1) + ETFox*QH2*v34_KcF/v34_Keq + ETFox*v34_KcF*v34_KmP2/v34_Keq + ETFred*Q*v34_KcR + ETFred*v34_KcR*v34_KmS2 + Q*QH2*v34_KcR*v34_KmS1/v34_KiP2 + Q*v34_KcR*v34_KmS1 + QH2*v34_KcF*v34_KmP1/v34_Keq, 2);
            break;
        case 206:
            dwdp[28] = -1.0*QH2*pow(v34_KcF, 2)*v34_KcR*v34_v34_ETF_QO*(-ETFox*QH2/v34_Keq + ETFred*Q)/(v34_Keq*pow(ETFox*ETFred*v34_KcF*v34_KmP2/(v34_Keq*v34_KiS1) + ETFox*QH2*v34_KcF/v34_Keq + ETFox*v34_KcF*v34_KmP2/v34_Keq + ETFred*Q*v34_KcR + ETFred*v34_KcR*v34_KmS2 + Q*QH2*v34_KcR*v34_KmS1/v34_KiP2 + Q*v34_KcR*v34_KmS1 + QH2*v34_KcF*v34_KmP1/v34_Keq, 2));
            break;
        case 207:
            dwdp[28] = -1.0*ETFred*v34_KcF*pow(v34_KcR, 2)*v34_v34_ETF_QO*(-ETFox*QH2/v34_Keq + ETFred*Q)/pow(ETFox*ETFred*v34_KcF*v34_KmP2/(v34_Keq*v34_KiS1) + ETFox*QH2*v34_KcF/v34_Keq + ETFox*v34_KcF*v34_KmP2/v34_Keq + ETFred*Q*v34_KcR + ETFred*v34_KcR*v34_KmS2 + Q*QH2*v34_KcR*v34_KmS1/v34_KiP2 + Q*v34_KcR*v34_KmS1 + QH2*v34_KcF*v34_KmP1/v34_Keq, 2);
            break;
        case 208:
            dwdp[28] = 1.0*v34_KcF*v34_KcR*v34_v34_ETF_QO*(-ETFox*QH2/v34_Keq + ETFred*Q)*(-Q*QH2*v34_KcR/v34_KiP2 - Q*v34_KcR)/pow(ETFox*ETFred*v34_KcF*v34_KmP2/(v34_Keq*v34_KiS1) + ETFox*QH2*v34_KcF/v34_Keq + ETFox*v34_KcF*v34_KmP2/v34_Keq + ETFred*Q*v34_KcR + ETFred*v34_KcR*v34_KmS2 + Q*QH2*v34_KcR*v34_KmS1/v34_KiP2 + Q*v34_KcR*v34_KmS1 + QH2*v34_KcF*v34_KmP1/v34_Keq, 2);
            break;
        case 209:
            dwdp[29] = 1.0*v35_KcF*v35_KcR*(ETFox*FADH2 - ETFred*FAD/v35_Keq)/(ETFox*ETFred*FAD*v35_KcF/(v35_Keq*v35_KiS2) + ETFox*ETFred*FADH2*v35_KcR/v35_KiP1 + ETFox*FAD*v35_KcR*v35_KmS1/v35_KiP2 + ETFox*FADH2*v35_KcR + ETFox*v35_KcR*v35_KmS1 + ETFred*FAD*v35_KcF/v35_Keq + ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/v35_Keq + FAD*v35_KcF*v35_KmP1/v35_Keq + FADH2*v35_KcR*v35_KmS2 + v35_KcR*v35_KiS1*v35_KmS2);
            break;
        case 210:
            dwdp[29] = 1.0*v35_KcF*v35_KcR*v35_v35_ACD*(ETFox*FADH2 - ETFred*FAD/v35_Keq)*(-ETFox*ETFred*FADH2/v35_KiP1 - ETFox*FAD*v35_KmS1/v35_KiP2 - ETFox*FADH2 - ETFox*v35_KmS1 - FADH2*v35_KmS2 - v35_KiS1*v35_KmS2)/pow(ETFox*ETFred*FAD*v35_KcF/(v35_Keq*v35_KiS2) + ETFox*ETFred*FADH2*v35_KcR/v35_KiP1 + ETFox*FAD*v35_KcR*v35_KmS1/v35_KiP2 + ETFox*FADH2*v35_KcR + ETFox*v35_KcR*v35_KmS1 + ETFred*FAD*v35_KcF/v35_Keq + ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/v35_Keq + FAD*v35_KcF*v35_KmP1/v35_Keq + FADH2*v35_KcR*v35_KmS2 + v35_KcR*v35_KiS1*v35_KmS2, 2) + 1.0*v35_KcF*v35_v35_ACD*(ETFox*FADH2 - ETFred*FAD/v35_Keq)/(ETFox*ETFred*FAD*v35_KcF/(v35_Keq*v35_KiS2) + ETFox*ETFred*FADH2*v35_KcR/v35_KiP1 + ETFox*FAD*v35_KcR*v35_KmS1/v35_KiP2 + ETFox*FADH2*v35_KcR + ETFox*v35_KcR*v35_KmS1 + ETFred*FAD*v35_KcF/v35_Keq + ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/v35_Keq + FAD*v35_KcF*v35_KmP1/v35_Keq + FADH2*v35_KcR*v35_KmS2 + v35_KcR*v35_KiS1*v35_KmS2);
            break;
        case 211:
            dwdp[29] = 1.0*v35_KcF*v35_KcR*v35_v35_ACD*(ETFox*FADH2 - ETFred*FAD/v35_Keq)*(-ETFox*ETFred*FAD/(v35_Keq*v35_KiS2) - ETFred*FAD/v35_Keq - ETFred*FADH2*v35_KmP2/(v35_Keq*v35_KiS1) - ETFred*v35_KmP2/v35_Keq - FAD*v35_KmP1/v35_Keq)/pow(ETFox*ETFred*FAD*v35_KcF/(v35_Keq*v35_KiS2) + ETFox*ETFred*FADH2*v35_KcR/v35_KiP1 + ETFox*FAD*v35_KcR*v35_KmS1/v35_KiP2 + ETFox*FADH2*v35_KcR + ETFox*v35_KcR*v35_KmS1 + ETFred*FAD*v35_KcF/v35_Keq + ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/v35_Keq + FAD*v35_KcF*v35_KmP1/v35_Keq + FADH2*v35_KcR*v35_KmS2 + v35_KcR*v35_KiS1*v35_KmS2, 2) + 1.0*v35_KcR*v35_v35_ACD*(ETFox*FADH2 - ETFred*FAD/v35_Keq)/(ETFox*ETFred*FAD*v35_KcF/(v35_Keq*v35_KiS2) + ETFox*ETFred*FADH2*v35_KcR/v35_KiP1 + ETFox*FAD*v35_KcR*v35_KmS1/v35_KiP2 + ETFox*FADH2*v35_KcR + ETFox*v35_KcR*v35_KmS1 + ETFred*FAD*v35_KcF/v35_Keq + ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/v35_Keq + FAD*v35_KcF*v35_KmP1/v35_Keq + FADH2*v35_KcR*v35_KmS2 + v35_KcR*v35_KiS1*v35_KmS2);
            break;
        case 212:
            dwdp[29] = 1.0*ETFred*FAD*v35_KcF*v35_KcR*v35_v35_ACD/(pow(v35_Keq, 2)*(ETFox*ETFred*FAD*v35_KcF/(v35_Keq*v35_KiS2) + ETFox*ETFred*FADH2*v35_KcR/v35_KiP1 + ETFox*FAD*v35_KcR*v35_KmS1/v35_KiP2 + ETFox*FADH2*v35_KcR + ETFox*v35_KcR*v35_KmS1 + ETFred*FAD*v35_KcF/v35_Keq + ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/v35_Keq + FAD*v35_KcF*v35_KmP1/v35_Keq + FADH2*v35_KcR*v35_KmS2 + v35_KcR*v35_KiS1*v35_KmS2)) + 1.0*v35_KcF*v35_KcR*v35_v35_ACD*(ETFox*FADH2 - ETFred*FAD/v35_Keq)*(ETFox*ETFred*FAD*v35_KcF/(pow(v35_Keq, 2)*v35_KiS2) + ETFred*FAD*v35_KcF/pow(v35_Keq, 2) + ETFred*FADH2*v35_KcF*v35_KmP2/(pow(v35_Keq, 2)*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/pow(v35_Keq, 2) + FAD*v35_KcF*v35_KmP1/pow(v35_Keq, 2))/pow(ETFox*ETFred*FAD*v35_KcF/(v35_Keq*v35_KiS2) + ETFox*ETFred*FADH2*v35_KcR/v35_KiP1 + ETFox*FAD*v35_KcR*v35_KmS1/v35_KiP2 + ETFox*FADH2*v35_KcR + ETFox*v35_KcR*v35_KmS1 + ETFred*FAD*v35_KcF/v35_Keq + ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/v35_Keq + FAD*v35_KcF*v35_KmP1/v35_Keq + FADH2*v35_KcR*v35_KmS2 + v35_KcR*v35_KiS1*v35_KmS2, 2);
            break;
        case 213:
            dwdp[29] = 1.0*ETFox*FAD*v35_KcF*pow(v35_KcR, 2)*v35_KmS1*v35_v35_ACD*(ETFox*FADH2 - ETFred*FAD/v35_Keq)/(pow(v35_KiP2, 2)*pow(ETFox*ETFred*FAD*v35_KcF/(v35_Keq*v35_KiS2) + ETFox*ETFred*FADH2*v35_KcR/v35_KiP1 + ETFox*FAD*v35_KcR*v35_KmS1/v35_KiP2 + ETFox*FADH2*v35_KcR + ETFox*v35_KcR*v35_KmS1 + ETFred*FAD*v35_KcF/v35_Keq + ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/v35_Keq + FAD*v35_KcF*v35_KmP1/v35_Keq + FADH2*v35_KcR*v35_KmS2 + v35_KcR*v35_KiS1*v35_KmS2, 2));
            break;
        case 214:
            dwdp[29] = 1.0*ETFox*ETFred*FADH2*v35_KcF*pow(v35_KcR, 2)*v35_v35_ACD*(ETFox*FADH2 - ETFred*FAD/v35_Keq)/(pow(v35_KiP1, 2)*pow(ETFox*ETFred*FAD*v35_KcF/(v35_Keq*v35_KiS2) + ETFox*ETFred*FADH2*v35_KcR/v35_KiP1 + ETFox*FAD*v35_KcR*v35_KmS1/v35_KiP2 + ETFox*FADH2*v35_KcR + ETFox*v35_KcR*v35_KmS1 + ETFred*FAD*v35_KcF/v35_Keq + ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/v35_Keq + FAD*v35_KcF*v35_KmP1/v35_Keq + FADH2*v35_KcR*v35_KmS2 + v35_KcR*v35_KiS1*v35_KmS2, 2));
            break;
        case 215:
            dwdp[29] = 1.0*ETFox*ETFred*FAD*pow(v35_KcF, 2)*v35_KcR*v35_v35_ACD*(ETFox*FADH2 - ETFred*FAD/v35_Keq)/(v35_Keq*pow(v35_KiS2, 2)*pow(ETFox*ETFred*FAD*v35_KcF/(v35_Keq*v35_KiS2) + ETFox*ETFred*FADH2*v35_KcR/v35_KiP1 + ETFox*FAD*v35_KcR*v35_KmS1/v35_KiP2 + ETFox*FADH2*v35_KcR + ETFox*v35_KcR*v35_KmS1 + ETFred*FAD*v35_KcF/v35_Keq + ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/v35_Keq + FAD*v35_KcF*v35_KmP1/v35_Keq + FADH2*v35_KcR*v35_KmS2 + v35_KcR*v35_KiS1*v35_KmS2, 2));
            break;
        case 216:
            dwdp[29] = 1.0*v35_KcF*v35_KcR*v35_v35_ACD*(ETFox*FADH2 - ETFred*FAD/v35_Keq)*(ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*pow(v35_KiS1, 2)) - v35_KcR*v35_KmS2)/pow(ETFox*ETFred*FAD*v35_KcF/(v35_Keq*v35_KiS2) + ETFox*ETFred*FADH2*v35_KcR/v35_KiP1 + ETFox*FAD*v35_KcR*v35_KmS1/v35_KiP2 + ETFox*FADH2*v35_KcR + ETFox*v35_KcR*v35_KmS1 + ETFred*FAD*v35_KcF/v35_Keq + ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/v35_Keq + FAD*v35_KcF*v35_KmP1/v35_Keq + FADH2*v35_KcR*v35_KmS2 + v35_KcR*v35_KiS1*v35_KmS2, 2);
            break;
        case 217:
            dwdp[29] = 1.0*v35_KcF*v35_KcR*v35_v35_ACD*(ETFox*FADH2 - ETFred*FAD/v35_Keq)*(-ETFred*FADH2*v35_KcF/(v35_Keq*v35_KiS1) - ETFred*v35_KcF/v35_Keq)/pow(ETFox*ETFred*FAD*v35_KcF/(v35_Keq*v35_KiS2) + ETFox*ETFred*FADH2*v35_KcR/v35_KiP1 + ETFox*FAD*v35_KcR*v35_KmS1/v35_KiP2 + ETFox*FADH2*v35_KcR + ETFox*v35_KcR*v35_KmS1 + ETFred*FAD*v35_KcF/v35_Keq + ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/v35_Keq + FAD*v35_KcF*v35_KmP1/v35_Keq + FADH2*v35_KcR*v35_KmS2 + v35_KcR*v35_KiS1*v35_KmS2, 2);
            break;
        case 218:
            dwdp[29] = -1.0*FAD*pow(v35_KcF, 2)*v35_KcR*v35_v35_ACD*(ETFox*FADH2 - ETFred*FAD/v35_Keq)/(v35_Keq*pow(ETFox*ETFred*FAD*v35_KcF/(v35_Keq*v35_KiS2) + ETFox*ETFred*FADH2*v35_KcR/v35_KiP1 + ETFox*FAD*v35_KcR*v35_KmS1/v35_KiP2 + ETFox*FADH2*v35_KcR + ETFox*v35_KcR*v35_KmS1 + ETFred*FAD*v35_KcF/v35_Keq + ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/v35_Keq + FAD*v35_KcF*v35_KmP1/v35_Keq + FADH2*v35_KcR*v35_KmS2 + v35_KcR*v35_KiS1*v35_KmS2, 2));
            break;
        case 219:
            dwdp[29] = 1.0*v35_KcF*v35_KcR*v35_v35_ACD*(ETFox*FADH2 - ETFred*FAD/v35_Keq)*(-FADH2*v35_KcR - v35_KcR*v35_KiS1)/pow(ETFox*ETFred*FAD*v35_KcF/(v35_Keq*v35_KiS2) + ETFox*ETFred*FADH2*v35_KcR/v35_KiP1 + ETFox*FAD*v35_KcR*v35_KmS1/v35_KiP2 + ETFox*FADH2*v35_KcR + ETFox*v35_KcR*v35_KmS1 + ETFred*FAD*v35_KcF/v35_Keq + ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/v35_Keq + FAD*v35_KcF*v35_KmP1/v35_Keq + FADH2*v35_KcR*v35_KmS2 + v35_KcR*v35_KiS1*v35_KmS2, 2);
            break;
        case 220:
            dwdp[29] = 1.0*v35_KcF*v35_KcR*v35_v35_ACD*(ETFox*FADH2 - ETFred*FAD/v35_Keq)*(-ETFox*FAD*v35_KcR/v35_KiP2 - ETFox*v35_KcR)/pow(ETFox*ETFred*FAD*v35_KcF/(v35_Keq*v35_KiS2) + ETFox*ETFred*FADH2*v35_KcR/v35_KiP1 + ETFox*FAD*v35_KcR*v35_KmS1/v35_KiP2 + ETFox*FADH2*v35_KcR + ETFox*v35_KcR*v35_KmS1 + ETFred*FAD*v35_KcF/v35_Keq + ETFred*FADH2*v35_KcF*v35_KmP2/(v35_Keq*v35_KiS1) + ETFred*v35_KcF*v35_KmP2/v35_Keq + FAD*v35_KcF*v35_KmP1/v35_Keq + FADH2*v35_KcR*v35_KmS2 + v35_KcR*v35_KiS1*v35_KmS2, 2);
            break;
        case 221:
            dwdp[30] = 1.0*v36_KcF*v36_KcR*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)/(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip);
            break;
        case 222:
            dwdp[30] = 1.0*v36_KcF*v36_KcR*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)*(-ADP*CO2*Pyr*v36_KmA/v36_Kiq - ADP*Pyr*v36_KmB - ATP*CO2*OXA*v36_KmA/v36_Kir - ATP*CO2*Pyr - ATP*CO2*v36_KmC - ATP*Pyr*v36_KmB - CO2*Pyr*v36_KmA - amici::pi*CO2*Pyr*v36_KmA/v36_Kip - Pyr*v36_Kia*v36_KmB - amici::pi*Pyr*v36_Kia*v36_KmB/v36_Kip)/pow(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip, 2) + 1.0*v36_KcF*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)/(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip);
            break;
        case 223:
            dwdp[30] = 1.0*v36_KcF*v36_KcR*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)*(-ADP*CO2*OXA*v36_KmP/(v36_Keq*v36_Kib) - ADP*OXA*v36_KmP/v36_Keq - amici::pi*ADP*OXA/v36_Keq - amici::pi*ADP*Pyr*v36_KmR/(v36_Keq*v36_Kic) - amici::pi*ADP*v36_KmR/v36_Keq - ATP*OXA*v36_KmP/v36_Keq - CO2*OXA*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) - amici::pi*CO2*OXA*v36_KmQ/(v36_Keq*v36_Kib) - OXA*v36_Kip*v36_KmQ/v36_Keq - amici::pi*OXA*v36_KmQ/v36_Keq)/pow(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip, 2) + 1.0*v36_KcR*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)/(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip);
            break;
        case 224:
            dwdp[30] = 1.0*ATP*CO2*OXA*v36_KcF*pow(v36_KcR, 2)*v36_KmA*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)/(pow(v36_Kir, 2)*pow(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip, 2));
            break;
        case 225:
            dwdp[30] = 1.0*v36_KcF*v36_KcR*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)*(ADP*CO2*Pyr*v36_KcR*v36_KmA/pow(v36_Kiq, 2) - CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib))/pow(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip, 2);
            break;
        case 226:
            dwdp[30] = 1.0*v36_KcF*v36_KcR*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)*(amici::pi*CO2*Pyr*v36_KcR*v36_KmA/pow(v36_Kip, 2) - OXA*v36_KcF*v36_KmQ/v36_Keq + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/pow(v36_Kip, 2))/pow(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip, 2);
            break;
        case 227:
            dwdp[30] = 1.0*amici::pi*ADP*Pyr*pow(v36_KcF, 2)*v36_KcR*v36_KmR*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)/(v36_Keq*pow(v36_Kic, 2)*pow(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip, 2));
            break;
        case 228:
            dwdp[30] = 1.0*v36_KcF*v36_KcR*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)*(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*pow(v36_Kib, 2)) + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*pow(v36_Kib, 2)) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*pow(v36_Kib, 2)))/pow(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip, 2);
            break;
        case 229:
            dwdp[30] = 1.0*v36_KcF*v36_KcR*v36_v36_PC*(-Pyr*v36_KcR*v36_KmB - amici::pi*Pyr*v36_KcR*v36_KmB/v36_Kip)*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)/pow(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip, 2);
            break;
        case 230:
            dwdp[30] = 1.0*amici::pi*ADP*OXA*v36_KcF*v36_KcR*v36_v36_PC/(pow(v36_Keq, 2)*(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip)) + 1.0*v36_KcF*v36_KcR*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)*(ADP*CO2*OXA*v36_KcF*v36_KmP/(pow(v36_Keq, 2)*v36_Kib) + ADP*OXA*v36_KcF*v36_KmP/pow(v36_Keq, 2) + amici::pi*ADP*OXA*v36_KcF/pow(v36_Keq, 2) + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(pow(v36_Keq, 2)*v36_Kic) + amici::pi*ADP*v36_KcF*v36_KmR/pow(v36_Keq, 2) + ATP*OXA*v36_KcF*v36_KmP/pow(v36_Keq, 2) + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(pow(v36_Keq, 2)*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(pow(v36_Keq, 2)*v36_Kib) + OXA*v36_KcF*v36_Kip*v36_KmQ/pow(v36_Keq, 2) + amici::pi*OXA*v36_KcF*v36_KmQ/pow(v36_Keq, 2))/pow(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip, 2);
            break;
        case 231:
            dwdp[30] = 1.0*v36_KcF*v36_KcR*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)*(-amici::pi*ADP*Pyr*v36_KcF/(v36_Keq*v36_Kic) - amici::pi*ADP*v36_KcF/v36_Keq)/pow(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip, 2);
            break;
        case 232:
            dwdp[30] = 1.0*v36_KcF*v36_KcR*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)*(-amici::pi*CO2*OXA*v36_KcF/(v36_Keq*v36_Kib) - OXA*v36_KcF*v36_Kip/v36_Keq - amici::pi*OXA*v36_KcF/v36_Keq)/pow(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip, 2);
            break;
        case 233:
            dwdp[30] = 1.0*v36_KcF*v36_KcR*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)*(-ADP*CO2*OXA*v36_KcF/(v36_Keq*v36_Kib) - ADP*OXA*v36_KcF/v36_Keq - ATP*OXA*v36_KcF/v36_Keq - CO2*OXA*v36_KcF*v36_Kiq/(v36_Keq*v36_Kib))/pow(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip, 2);
            break;
        case 234:
            dwdp[30] = -1.0*ATP*CO2*v36_KcF*pow(v36_KcR, 2)*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)/pow(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip, 2);
            break;
        case 235:
            dwdp[30] = 1.0*v36_KcF*v36_KcR*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)*(-ADP*Pyr*v36_KcR - ATP*Pyr*v36_KcR - Pyr*v36_KcR*v36_Kia - amici::pi*Pyr*v36_KcR*v36_Kia/v36_Kip)/pow(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip, 2);
            break;
        case 236:
            dwdp[30] = 1.0*v36_KcF*v36_KcR*v36_v36_PC*(-amici::pi*ADP*OXA/v36_Keq + ATP*CO2*Pyr)*(-ADP*CO2*Pyr*v36_KcR/v36_Kiq - ATP*CO2*OXA*v36_KcR/v36_Kir - CO2*Pyr*v36_KcR - amici::pi*CO2*Pyr*v36_KcR/v36_Kip)/pow(ADP*CO2*OXA*v36_KcF*v36_KmP/(v36_Keq*v36_Kib) + ADP*CO2*Pyr*v36_KcR*v36_KmA/v36_Kiq + ADP*OXA*v36_KcF*v36_KmP/v36_Keq + amici::pi*ADP*OXA*v36_KcF/v36_Keq + amici::pi*ADP*Pyr*v36_KcF*v36_KmR/(v36_Keq*v36_Kic) + ADP*Pyr*v36_KcR*v36_KmB + amici::pi*ADP*v36_KcF*v36_KmR/v36_Keq + ATP*CO2*OXA*v36_KcR*v36_KmA/v36_Kir + ATP*CO2*Pyr*v36_KcR + ATP*CO2*v36_KcR*v36_KmC + ATP*OXA*v36_KcF*v36_KmP/v36_Keq + ATP*Pyr*v36_KcR*v36_KmB + CO2*OXA*v36_KcF*v36_Kiq*v36_KmP/(v36_Keq*v36_Kib) + amici::pi*CO2*OXA*v36_KcF*v36_KmQ/(v36_Keq*v36_Kib) + CO2*Pyr*v36_KcR*v36_KmA + amici::pi*CO2*Pyr*v36_KcR*v36_KmA/v36_Kip + OXA*v36_KcF*v36_Kip*v36_KmQ/v36_Keq + amici::pi*OXA*v36_KcF*v36_KmQ/v36_Keq + Pyr*v36_KcR*v36_Kia*v36_KmB + amici::pi*Pyr*v36_KcR*v36_Kia*v36_KmB/v36_Kip, 2);
            break;
        case 237:
            dwdp[31] = 1.0*G3P*v37_V/(G3P + v37_K);
            break;
        case 238:
            dwdp[31] = 1.0*G3P*v37_v37_GUT2P/(G3P + v37_K);
            break;
        case 239:
            dwdp[31] = -1.0*G3P*v37_V*v37_v37_GUT2P/pow(G3P + v37_K, 2);
            break;
        case 240:
            dwdp[32] = 1.0*NADH_cyt*v38_V/(NADH_cyt + v38_K);
            break;
        case 241:
            dwdp[32] = 1.0*NADH_cyt*v38_v38_GUT2P/(NADH_cyt + v38_K);
            break;
        case 242:
            dwdp[32] = -1.0*NADH_cyt*v38_V*v38_v38_GUT2P/pow(NADH_cyt + v38_K, 2);
            break;
        case 243:
            dwdp[33] = 1.0*Mal_cyt*NADP_cyt*v39_Kcat/((Mal_cyt + v39_Kmal)*(NADP_cyt + v39_Knadp));
            break;
        case 244:
            dwdp[33] = -1.0*Mal_cyt*NADP_cyt*v39_Kcat*v39_v39_MDH/((Mal_cyt + v39_Kmal)*pow(NADP_cyt + v39_Knadp, 2));
            break;
        case 245:
            dwdp[33] = -1.0*Mal_cyt*NADP_cyt*v39_Kcat*v39_v39_MDH/(pow(Mal_cyt + v39_Kmal, 2)*(NADP_cyt + v39_Knadp));
            break;
        case 246:
            dwdp[33] = 1.0*Mal_cyt*NADP_cyt*v39_v39_MDH/((Mal_cyt + v39_Kmal)*(NADP_cyt + v39_Knadp));
            break;
        case 247:
            dwdp[34] = -1.0*GAP*NAD*v4_V4/((GAP + v4_K4GAP)*pow(NAD + v4_K4NAD, 2));
            break;
        case 248:
            dwdp[34] = -1.0*GAP*NAD*v4_V4/(pow(GAP + v4_K4GAP, 2)*(NAD + v4_K4NAD));
            break;
        case 249:
            dwdp[34] = 1.0*GAP*NAD/((GAP + v4_K4GAP)*(NAD + v4_K4NAD));
            break;
        case 250:
            dwdp[35] = 1.0*ADP_cyt*v40_V/(ADP_cyt + v40_K);
            break;
        case 251:
            dwdp[35] = -1.0*ADP_cyt*v40_V*v40_v40_AAC/pow(ADP_cyt + v40_K, 2);
            break;
        case 252:
            dwdp[35] = 1.0*ADP_cyt*v40_v40_AAC/(ADP_cyt + v40_K);
            break;
        case 253:
            dwdp[36] = -1.0*CO2*NADPH_cyt*OG_cyt/(CO2*NADPH_cyt*OG_cyt*v41_phir0 + CO2*NADPH_cyt*v41_phir1 + CO2*OG_cyt*v41_phir2 + CO2*v41_phir12 + NADPH_cyt*OG_cyt*v41_phir3 + NADPH_cyt*v41_phir13 + OG_cyt*v41_phir23 + v41_phir123) + 1.0*IsoCitcyt*NADP_cyt/(IsoCitcyt*NADP_cyt*v41_phi0 + IsoCitcyt*v41_phi2 + NADP_cyt*v41_phi1 + v41_phi12);
            break;
        case 254:
            dwdp[36] = 1.0*CO2*NADPH_cyt*OG_cyt*v41_v41_IDHc/pow(CO2*NADPH_cyt*OG_cyt*v41_phir0 + CO2*NADPH_cyt*v41_phir1 + CO2*OG_cyt*v41_phir2 + CO2*v41_phir12 + NADPH_cyt*OG_cyt*v41_phir3 + NADPH_cyt*v41_phir13 + OG_cyt*v41_phir23 + v41_phir123, 2);
            break;
        case 255:
            dwdp[36] = 1.0*CO2*NADPH_cyt*pow(OG_cyt, 2)*v41_v41_IDHc/pow(CO2*NADPH_cyt*OG_cyt*v41_phir0 + CO2*NADPH_cyt*v41_phir1 + CO2*OG_cyt*v41_phir2 + CO2*v41_phir12 + NADPH_cyt*OG_cyt*v41_phir3 + NADPH_cyt*v41_phir13 + OG_cyt*v41_phir23 + v41_phir123, 2);
            break;
        case 256:
            dwdp[36] = 1.0*CO2*pow(NADPH_cyt, 2)*OG_cyt*v41_v41_IDHc/pow(CO2*NADPH_cyt*OG_cyt*v41_phir0 + CO2*NADPH_cyt*v41_phir1 + CO2*OG_cyt*v41_phir2 + CO2*v41_phir12 + NADPH_cyt*OG_cyt*v41_phir3 + NADPH_cyt*v41_phir13 + OG_cyt*v41_phir23 + v41_phir123, 2);
            break;
        case 257:
            dwdp[36] = 1.0*pow(CO2, 2)*NADPH_cyt*OG_cyt*v41_v41_IDHc/pow(CO2*NADPH_cyt*OG_cyt*v41_phir0 + CO2*NADPH_cyt*v41_phir1 + CO2*OG_cyt*v41_phir2 + CO2*v41_phir12 + NADPH_cyt*OG_cyt*v41_phir3 + NADPH_cyt*v41_phir13 + OG_cyt*v41_phir23 + v41_phir123, 2);
            break;
        case 258:
            dwdp[36] = 1.0*CO2*pow(NADPH_cyt, 2)*pow(OG_cyt, 2)*v41_v41_IDHc/pow(CO2*NADPH_cyt*OG_cyt*v41_phir0 + CO2*NADPH_cyt*v41_phir1 + CO2*OG_cyt*v41_phir2 + CO2*v41_phir12 + NADPH_cyt*OG_cyt*v41_phir3 + NADPH_cyt*v41_phir13 + OG_cyt*v41_phir23 + v41_phir123, 2);
            break;
        case 259:
            dwdp[36] = 1.0*pow(CO2, 2)*NADPH_cyt*pow(OG_cyt, 2)*v41_v41_IDHc/pow(CO2*NADPH_cyt*OG_cyt*v41_phir0 + CO2*NADPH_cyt*v41_phir1 + CO2*OG_cyt*v41_phir2 + CO2*v41_phir12 + NADPH_cyt*OG_cyt*v41_phir3 + NADPH_cyt*v41_phir13 + OG_cyt*v41_phir23 + v41_phir123, 2);
            break;
        case 260:
            dwdp[36] = 1.0*pow(CO2, 2)*pow(NADPH_cyt, 2)*OG_cyt*v41_v41_IDHc/pow(CO2*NADPH_cyt*OG_cyt*v41_phir0 + CO2*NADPH_cyt*v41_phir1 + CO2*OG_cyt*v41_phir2 + CO2*v41_phir12 + NADPH_cyt*OG_cyt*v41_phir3 + NADPH_cyt*v41_phir13 + OG_cyt*v41_phir23 + v41_phir123, 2);
            break;
        case 261:
            dwdp[36] = 1.0*pow(CO2, 2)*pow(NADPH_cyt, 2)*pow(OG_cyt, 2)*v41_v41_IDHc/pow(CO2*NADPH_cyt*OG_cyt*v41_phir0 + CO2*NADPH_cyt*v41_phir1 + CO2*OG_cyt*v41_phir2 + CO2*v41_phir12 + NADPH_cyt*OG_cyt*v41_phir3 + NADPH_cyt*v41_phir13 + OG_cyt*v41_phir23 + v41_phir123, 2);
            break;
        case 262:
            dwdp[36] = -1.0*IsoCitcyt*NADP_cyt*v41_v41_IDHc/pow(IsoCitcyt*NADP_cyt*v41_phi0 + IsoCitcyt*v41_phi2 + NADP_cyt*v41_phi1 + v41_phi12, 2);
            break;
        case 263:
            dwdp[36] = -1.0*pow(IsoCitcyt, 2)*NADP_cyt*v41_v41_IDHc/pow(IsoCitcyt*NADP_cyt*v41_phi0 + IsoCitcyt*v41_phi2 + NADP_cyt*v41_phi1 + v41_phi12, 2);
            break;
        case 264:
            dwdp[36] = -1.0*IsoCitcyt*pow(NADP_cyt, 2)*v41_v41_IDHc/pow(IsoCitcyt*NADP_cyt*v41_phi0 + IsoCitcyt*v41_phi2 + NADP_cyt*v41_phi1 + v41_phi12, 2);
            break;
        case 265:
            dwdp[36] = -1.0*pow(IsoCitcyt, 2)*pow(NADP_cyt, 2)*v41_v41_IDHc/pow(IsoCitcyt*NADP_cyt*v41_phi0 + IsoCitcyt*v41_phi2 + NADP_cyt*v41_phi1 + v41_phi12, 2);
            break;
        case 266:
            dwdp[37] = 1.0*(-IsoCit*Mal_cyt*v42_KcR/(v42_KiP1*v42_KiP2*v42_beta) + IsoCitcyt*Mal*v42_KcF/(v42_KiS1*v42_KiS2*v42_alpha))/(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1);
            break;
        case 267:
            dwdp[37] = 1.0*IsoCitcyt*Mal_cyt*v42_v42_CIC*(-IsoCit*Mal_cyt*v42_KcR/(v42_KiP1*v42_KiP2*v42_beta) + IsoCitcyt*Mal*v42_KcF/(v42_KiS1*v42_KiS2*v42_alpha))/(v42_KiP1*v42_KiS1*pow(v42_delta, 2)*pow(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1, 2));
            break;
        case 268:
            dwdp[37] = 1.0*IsoCit*Mal*v42_v42_CIC*(-IsoCit*Mal_cyt*v42_KcR/(v42_KiP1*v42_KiP2*v42_beta) + IsoCitcyt*Mal*v42_KcF/(v42_KiS1*v42_KiS2*v42_alpha))/(v42_KiP2*v42_KiS2*pow(v42_gamma, 2)*pow(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1, 2));
            break;
        case 269:
            dwdp[37] = 1.0*IsoCit*Mal_cyt*v42_KcR*v42_v42_CIC/(v42_KiP1*v42_KiP2*pow(v42_beta, 2)*(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1)) + 1.0*IsoCit*Mal_cyt*v42_v42_CIC*(-IsoCit*Mal_cyt*v42_KcR/(v42_KiP1*v42_KiP2*v42_beta) + IsoCitcyt*Mal*v42_KcF/(v42_KiS1*v42_KiS2*v42_alpha))/(v42_KiP1*v42_KiP2*pow(v42_beta, 2)*pow(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1, 2));
            break;
        case 270:
            dwdp[37] = -1.0*IsoCitcyt*Mal*v42_KcF*v42_v42_CIC/(v42_KiS1*v42_KiS2*pow(v42_alpha, 2)*(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1)) + 1.0*IsoCitcyt*Mal*v42_v42_CIC*(-IsoCit*Mal_cyt*v42_KcR/(v42_KiP1*v42_KiP2*v42_beta) + IsoCitcyt*Mal*v42_KcF/(v42_KiS1*v42_KiS2*v42_alpha))/(v42_KiS1*v42_KiS2*pow(v42_alpha, 2)*pow(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1, 2));
            break;
        case 271:
            dwdp[37] = -1.0*IsoCit*Mal_cyt*v42_v42_CIC/(v42_KiP1*v42_KiP2*v42_beta*(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1));
            break;
        case 272:
            dwdp[37] = 1.0*IsoCitcyt*Mal*v42_v42_CIC/(v42_KiS1*v42_KiS2*v42_alpha*(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1));
            break;
        case 273:
            dwdp[37] = 1.0*IsoCit*Mal_cyt*v42_KcR*v42_v42_CIC/(v42_KiP1*pow(v42_KiP2, 2)*v42_beta*(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1)) + 1.0*v42_v42_CIC*(-IsoCit*Mal_cyt*v42_KcR/(v42_KiP1*v42_KiP2*v42_beta) + IsoCitcyt*Mal*v42_KcF/(v42_KiS1*v42_KiS2*v42_alpha))*(IsoCit*Mal/(pow(v42_KiP2, 2)*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*pow(v42_KiP2, 2)*v42_beta) + IsoCit/pow(v42_KiP2, 2))/pow(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1, 2);
            break;
        case 274:
            dwdp[37] = 1.0*IsoCit*Mal_cyt*v42_KcR*v42_v42_CIC/(pow(v42_KiP1, 2)*v42_KiP2*v42_beta*(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1)) + 1.0*v42_v42_CIC*(-IsoCit*Mal_cyt*v42_KcR/(v42_KiP1*v42_KiP2*v42_beta) + IsoCitcyt*Mal*v42_KcF/(v42_KiS1*v42_KiS2*v42_alpha))*(IsoCit*Mal_cyt/(pow(v42_KiP1, 2)*v42_KiP2*v42_beta) + IsoCitcyt*Mal_cyt/(pow(v42_KiP1, 2)*v42_KiS1*v42_delta) + Mal_cyt/pow(v42_KiP1, 2))/pow(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1, 2);
            break;
        case 275:
            dwdp[37] = -1.0*IsoCitcyt*Mal*v42_KcF*v42_v42_CIC/(v42_KiS1*pow(v42_KiS2, 2)*v42_alpha*(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1)) + 1.0*v42_v42_CIC*(-IsoCit*Mal_cyt*v42_KcR/(v42_KiP1*v42_KiP2*v42_beta) + IsoCitcyt*Mal*v42_KcF/(v42_KiS1*v42_KiS2*v42_alpha))*(IsoCit*Mal/(v42_KiP2*pow(v42_KiS2, 2)*v42_gamma) + IsoCitcyt*Mal/(v42_KiS1*pow(v42_KiS2, 2)*v42_alpha) + Mal/pow(v42_KiS2, 2))/pow(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1, 2);
            break;
        case 276:
            dwdp[37] = -1.0*IsoCitcyt*Mal*v42_KcF*v42_v42_CIC/(pow(v42_KiS1, 2)*v42_KiS2*v42_alpha*(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1)) + 1.0*v42_v42_CIC*(-IsoCit*Mal_cyt*v42_KcR/(v42_KiP1*v42_KiP2*v42_beta) + IsoCitcyt*Mal*v42_KcF/(v42_KiS1*v42_KiS2*v42_alpha))*(IsoCitcyt*Mal/(pow(v42_KiS1, 2)*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*pow(v42_KiS1, 2)*v42_delta) + IsoCitcyt/pow(v42_KiS1, 2))/pow(IsoCit*Mal/(v42_KiP2*v42_KiS2*v42_gamma) + IsoCit*Mal_cyt/(v42_KiP1*v42_KiP2*v42_beta) + IsoCit/v42_KiP2 + IsoCitcyt*Mal/(v42_KiS1*v42_KiS2*v42_alpha) + IsoCitcyt*Mal_cyt/(v42_KiP1*v42_KiS1*v42_delta) + IsoCitcyt/v42_KiS1 + Mal/v42_KiS2 + Mal_cyt/v42_KiP1 + 1, 2);
            break;
        case 277:
            dwdp[38] = 1.0*ATP*v43_V/(ATP + v43_K);
            break;
        case 278:
            dwdp[38] = -1.0*ATP*v43_V*v43_v43_AAC/pow(ATP + v43_K, 2);
            break;
        case 279:
            dwdp[38] = 1.0*ATP*v43_v43_AAC/(ATP + v43_K);
            break;
        case 280:
            dwdp[39] = 1.0*Mal*v44_Kcat/(Mal + v44_Km);
            break;
        case 281:
            dwdp[39] = -1.0*Mal*v44_Kcat*v44_v44_MDH/pow(Mal + v44_Km, 2);
            break;
        case 282:
            dwdp[39] = 1.0*Mal*v44_v44_MDH/(Mal + v44_Km);
            break;
        case 283:
            dwdp[40] = -1.0*ATP_cyt*PEP;
            break;
        case 284:
            dwdp[40] = 1.0*ADP_cyt*DPG;
            break;
        case 285:
            dwdp[41] = -1.0*ADP_cyt*PEP*v6_V6/(pow(ADP_cyt + v6_K6ADP, 2)*(PEP + v6_K6PEP));
            break;
        case 286:
            dwdp[41] = -1.0*ADP_cyt*PEP*v6_V6/((ADP_cyt + v6_K6ADP)*pow(PEP + v6_K6PEP, 2));
            break;
        case 287:
            dwdp[41] = 1.0*ADP_cyt*PEP/((ADP_cyt + v6_K6ADP)*(PEP + v6_K6PEP));
            break;
        case 288:
            dwdp[42] = -1.0*LAC*NAD;
            break;
        case 289:
            dwdp[42] = 1.0*NADH_cyt*PYR_cyt;
            break;
        case 290:
            dwdp[43] = 1.0*PYR_cyt*v8_V/(PYR_cyt + v8_K);
            break;
        case 291:
            dwdp[43] = -1.0*PYR_cyt*v8_V*v8_v8_PYC/pow(PYR_cyt + v8_K, 2);
            break;
        case 292:
            dwdp[43] = 1.0*PYR_cyt*v8_v8_PYC/(PYR_cyt + v8_K);
            break;
        case 293:
            dwdp[44] = 1.0*CoA*NAD_p*Pyr*v9_KcF/(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NAD_p*Pyr*v9_KmB/v9_Kiq + CoA*NADH*Pyr*v9_KmC/v9_Kir + CoA*NAD_p*Pyr + CoA*NAD_p*v9_KmA + CoA*Pyr*v9_KmC + NAD_p*Pyr*v9_KmB);
            break;
        case 294:
            dwdp[44] = 1.0*CoA*NAD_p*Pyr*v9_v9_PDC/(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NAD_p*Pyr*v9_KmB/v9_Kiq + CoA*NADH*Pyr*v9_KmC/v9_Kir + CoA*NAD_p*Pyr + CoA*NAD_p*v9_KmA + CoA*Pyr*v9_KmC + NAD_p*Pyr*v9_KmB);
            break;
        case 295:
            dwdp[44] = 1.0*pow(CoA, 2)*NADH*NAD_p*pow(Pyr, 2)*v9_KcF*v9_KmC*v9_v9_PDC/(pow(v9_Kir, 2)*pow(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NAD_p*Pyr*v9_KmB/v9_Kiq + CoA*NADH*Pyr*v9_KmC/v9_Kir + CoA*NAD_p*Pyr + CoA*NAD_p*v9_KmA + CoA*Pyr*v9_KmC + NAD_p*Pyr*v9_KmB, 2));
            break;
        case 296:
            dwdp[44] = 1.0*CoA*NAD_p*Pyr*v9_KcF*v9_v9_PDC*(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*pow(v9_Kiq, 2)*v9_KmR) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*pow(v9_Kiq, 2)*v9_KmR) + Acetyl_CoA*NAD_p*Pyr*v9_KmB/pow(v9_Kiq, 2))/pow(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NAD_p*Pyr*v9_KmB/v9_Kiq + CoA*NADH*Pyr*v9_KmC/v9_Kir + CoA*NAD_p*Pyr + CoA*NAD_p*v9_KmA + CoA*Pyr*v9_KmC + NAD_p*Pyr*v9_KmB, 2);
            break;
        case 297:
            dwdp[44] = 1.0*CoA*NAD_p*Pyr*v9_KcF*v9_v9_PDC*(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*pow(v9_Kip, 2)*v9_Kiq*v9_KmR) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(pow(v9_Kip, 2)*v9_Kiq*v9_KmR))/pow(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NAD_p*Pyr*v9_KmB/v9_Kiq + CoA*NADH*Pyr*v9_KmC/v9_Kir + CoA*NAD_p*Pyr + CoA*NAD_p*v9_KmA + CoA*Pyr*v9_KmC + NAD_p*Pyr*v9_KmB, 2);
            break;
        case 298:
            dwdp[44] = 1.0*CoA*NAD_p*Pyr*v9_KcF*v9_v9_PDC*(-Acetyl_CoA*NADH*Pyr*v9_Kib*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) - Acetyl_CoA*NADH*v9_Kib*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR))/pow(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NAD_p*Pyr*v9_KmB/v9_Kiq + CoA*NADH*Pyr*v9_KmC/v9_Kir + CoA*NAD_p*Pyr + CoA*NAD_p*v9_KmA + CoA*Pyr*v9_KmC + NAD_p*Pyr*v9_KmB, 2);
            break;
        case 299:
            dwdp[44] = 1.0*CoA*NAD_p*Pyr*v9_KcF*v9_v9_PDC*(-Acetyl_CoA*NADH*Pyr*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) - Acetyl_CoA*NADH*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR))/pow(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NAD_p*Pyr*v9_KmB/v9_Kiq + CoA*NADH*Pyr*v9_KmC/v9_Kir + CoA*NAD_p*Pyr + CoA*NAD_p*v9_KmA + CoA*Pyr*v9_KmC + NAD_p*Pyr*v9_KmB, 2);
            break;
        case 300:
            dwdp[44] = 1.0*Acetyl_CoA*CoA*NADH*NAD_p*pow(Pyr, 2)*v9_KcF*v9_Kib*v9_Kic*v9_KmA*v9_KmP*v9_v9_PDC/(pow(v9_Kia, 2)*v9_Kip*v9_Kiq*v9_KmR*pow(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NAD_p*Pyr*v9_KmB/v9_Kiq + CoA*NADH*Pyr*v9_KmC/v9_Kir + CoA*NAD_p*Pyr + CoA*NAD_p*v9_KmA + CoA*Pyr*v9_KmC + NAD_p*Pyr*v9_KmB, 2));
            break;
        case 301:
            dwdp[44] = 1.0*CoA*NAD_p*Pyr*v9_KcF*v9_v9_PDC*(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*pow(v9_KmR, 2)) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*pow(v9_KmR, 2)))/pow(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NAD_p*Pyr*v9_KmB/v9_Kiq + CoA*NADH*Pyr*v9_KmC/v9_Kir + CoA*NAD_p*Pyr + CoA*NAD_p*v9_KmA + CoA*Pyr*v9_KmC + NAD_p*Pyr*v9_KmB, 2);
            break;
        case 302:
            dwdp[44] = 1.0*CoA*NAD_p*Pyr*v9_KcF*v9_v9_PDC*(-Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) - Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA/(v9_Kip*v9_Kiq*v9_KmR))/pow(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NAD_p*Pyr*v9_KmB/v9_Kiq + CoA*NADH*Pyr*v9_KmC/v9_Kir + CoA*NAD_p*Pyr + CoA*NAD_p*v9_KmA + CoA*Pyr*v9_KmC + NAD_p*Pyr*v9_KmB, 2);
            break;
        case 303:
            dwdp[44] = 1.0*CoA*NAD_p*Pyr*v9_KcF*v9_v9_PDC*(-CoA*NADH*Pyr/v9_Kir - CoA*Pyr)/pow(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NAD_p*Pyr*v9_KmB/v9_Kiq + CoA*NADH*Pyr*v9_KmC/v9_Kir + CoA*NAD_p*Pyr + CoA*NAD_p*v9_KmA + CoA*Pyr*v9_KmC + NAD_p*Pyr*v9_KmB, 2);
            break;
        case 304:
            dwdp[44] = 1.0*CoA*NAD_p*Pyr*v9_KcF*v9_v9_PDC*(-Acetyl_CoA*NAD_p*Pyr/v9_Kiq - NAD_p*Pyr)/pow(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NAD_p*Pyr*v9_KmB/v9_Kiq + CoA*NADH*Pyr*v9_KmC/v9_Kir + CoA*NAD_p*Pyr + CoA*NAD_p*v9_KmA + CoA*Pyr*v9_KmC + NAD_p*Pyr*v9_KmB, 2);
            break;
        case 305:
            dwdp[44] = 1.0*CoA*NAD_p*Pyr*v9_KcF*v9_v9_PDC*(-Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) - Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR) - CoA*NAD_p)/pow(Acetyl_CoA*NADH*Pyr*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kia*v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NADH*v9_Kib*v9_Kic*v9_KmA*v9_KmP/(v9_Kip*v9_Kiq*v9_KmR) + Acetyl_CoA*NAD_p*Pyr*v9_KmB/v9_Kiq + CoA*NADH*Pyr*v9_KmC/v9_Kir + CoA*NAD_p*Pyr + CoA*NAD_p*v9_KmA + CoA*Pyr*v9_KmC + NAD_p*Pyr*v9_KmB, 2);
            break;
    }
}