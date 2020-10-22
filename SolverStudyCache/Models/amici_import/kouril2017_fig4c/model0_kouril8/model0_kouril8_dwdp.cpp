#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_kouril8(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[3] = BPG*protGAPdh*(BPG*Vmfor*nadph/(KBPG*Knadph) - Vmarev*gap*nadp*pi/(Kgap*Knadp*Kpi))/(pow(KBPG, 2)*pow(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi), 2)*(1 + nadph/Knadph + nadp/Knadp)) - BPG*Vmfor*nadph*protGAPdh/(pow(KBPG, 2)*Knadph*(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*(1 + nadph/Knadph + nadp/Knadp));
            break;
        case 1:
            dwdp[3] = gap*protGAPdh*(1 + pi/Kpi)*(BPG*Vmfor*nadph/(KBPG*Knadph) - Vmarev*gap*nadp*pi/(Kgap*Knadp*Kpi))/(pow(Kgap, 2)*pow(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi), 2)*(1 + nadph/Knadph + nadp/Knadp)) + Vmarev*gap*nadp*pi*protGAPdh/(pow(Kgap, 2)*Knadp*Kpi*(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*(1 + nadph/Knadph + nadp/Knadp));
            break;
        case 2:
            dwdp[0] = ADP*protPGK*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(pow(KiADP, 2)*pow(ADP/KiADP + 1, 2)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G));
            break;
        case 3:
            dwdp[3] = nadp*protGAPdh*(BPG*Vmfor*nadph/(KBPG*Knadph) - Vmarev*gap*nadp*pi/(Kgap*Knadp*Kpi))/(pow(Knadp, 2)*(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*pow(1 + nadph/Knadph + nadp/Knadp, 2)) + Vmarev*gap*nadp*pi*protGAPdh/(Kgap*pow(Knadp, 2)*Kpi*(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*(1 + nadph/Knadph + nadp/Knadp));
            break;
        case 4:
            dwdp[3] = -BPG*Vmfor*nadph*protGAPdh/(KBPG*pow(Knadph, 2)*(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*(1 + nadph/Knadph + nadp/Knadp)) + nadph*protGAPdh*(BPG*Vmfor*nadph/(KBPG*Knadph) - Vmarev*gap*nadp*pi/(Kgap*Knadp*Kpi))/(pow(Knadph, 2)*(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*pow(1 + nadph/Knadph + nadp/Knadp, 2));
            break;
        case 5:
            dwdp[0] = ADP*BPG*VmrPGK*protPGK/(pow(KpgkADP, 2)*KpgkBPG*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G)) + ADP*BPG*protPGK*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(pow(KpgkADP, 2)*KpgkBPG*(ADP/KiADP + 1)*pow(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G, 2));
            break;
        case 6:
            dwdp[0] = -ATP*P3G*VmfPGK*protPGK/(pow(KpgkATP, 2)*KpgkP3G*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G)) + ATP*P3G*protPGK*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(pow(KpgkATP, 2)*KpgkP3G*(ADP/KiADP + 1)*pow(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G, 2));
            break;
        case 7:
            dwdp[0] = ADP*BPG*VmrPGK*protPGK/(KpgkADP*pow(KpgkBPG, 2)*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G)) + BPG*protPGK*(ADP/KpgkADP + 1)*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(pow(KpgkBPG, 2)*(ADP/KiADP + 1)*pow(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G, 2));
            break;
        case 8:
            dwdp[0] = -ATP*P3G*VmfPGK*protPGK/(KpgkATP*pow(KpgkP3G, 2)*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G)) + P3G*protPGK*(ATP/KpgkATP + 1)*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(pow(KpgkP3G, 2)*(ADP/KiADP + 1)*pow(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G, 2));
            break;
        case 9:
            dwdp[3] = pi*protGAPdh*(1 + gap/Kgap)*(BPG*Vmfor*nadph/(KBPG*Knadph) - Vmarev*gap*nadp*pi/(Kgap*Knadp*Kpi))/(pow(Kpi, 2)*pow(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi), 2)*(1 + nadph/Knadph + nadp/Knadp)) + Vmarev*gap*nadp*pi*protGAPdh/(Kgap*Knadp*pow(Kpi, 2)*(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*(1 + nadph/Knadph + nadp/Knadp));
            break;
        case 10:
            dwdp[3] = -gap*nadp*pi*protGAPdh/(Kgap*Knadp*Kpi*(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*(1 + nadph/Knadph + nadp/Knadp));
            break;
        case 11:
            dwdp[0] = ATP*P3G*protPGK/(KpgkATP*KpgkP3G*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G));
            break;
        case 12:
            dwdp[3] = BPG*nadph*protGAPdh/(KBPG*Knadph*(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*(1 + nadph/Knadph + nadp/Knadp));
            break;
        case 13:
            dwdp[0] = -ADP*BPG*protPGK/(KpgkADP*KpgkBPG*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G));
            break;
        case 14:
            dwdp[4] = glc*nadp;
            break;
        case 15:
            dwdp[1] = ADP*pep;
            break;
        case 16:
            dwdp[2] = BPG;
            break;
        case 17:
            dwdp[5] = gap;
            break;
        case 18:
            dwdp[3] = (BPG*Vmfor*nadph/(KBPG*Knadph) - Vmarev*gap*nadp*pi/(Kgap*Knadp*Kpi))/((BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*(1 + nadph/Knadph + nadp/Knadp));
            break;
        case 19:
            dwdp[0] = (-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/((ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G));
            break;
    }
}