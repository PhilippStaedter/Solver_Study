#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_kouril8(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -BPG*VmrPGK*protPGK/(KpgkADP*KpgkBPG*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G)) - BPG*protPGK*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(KpgkADP*KpgkBPG*(ADP/KiADP + 1)*pow(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G, 2)) - protPGK*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(KiADP*pow(ADP/KiADP + 1, 2)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G));
    dwdx[1] = kPK*pep;
    dwdx[2] = P3G*VmfPGK*protPGK/(KpgkATP*KpgkP3G*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G)) - P3G*protPGK*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(KpgkATP*KpgkP3G*(ADP/KiADP + 1)*pow(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G, 2));
    dwdx[3] = -ADP*VmrPGK*protPGK/(KpgkADP*KpgkBPG*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G)) - protPGK*(ADP/KpgkADP + 1)*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(KpgkBPG*(ADP/KiADP + 1)*pow(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G, 2));
    dwdx[4] = kdbpg;
    dwdx[5] = -protGAPdh*(BPG*Vmfor*nadph/(KBPG*Knadph) - Vmarev*gap*nadp*pi/(Kgap*Knadp*Kpi))/(KBPG*pow(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi), 2)*(1 + nadph/Knadph + nadp/Knadp)) + Vmfor*nadph*protGAPdh/(KBPG*Knadph*(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*(1 + nadph/Knadph + nadp/Knadp));
    dwdx[6] = ATP*VmfPGK*protPGK/(KpgkATP*KpgkP3G*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G)) - protPGK*(ATP/KpgkATP + 1)*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(KpgkP3G*(ADP/KiADP + 1)*pow(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G, 2));
    dwdx[7] = -protGAPdh*(1 + pi/Kpi)*(BPG*Vmfor*nadph/(KBPG*Knadph) - Vmarev*gap*nadp*pi/(Kgap*Knadp*Kpi))/(Kgap*pow(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi), 2)*(1 + nadph/Knadph + nadp/Knadp)) - Vmarev*nadp*pi*protGAPdh/(Kgap*Knadp*Kpi*(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*(1 + nadph/Knadph + nadp/Knadp));
    dwdx[8] = kdgap;
    dwdx[9] = kGDH*nadp;
    dwdx[10] = -protGAPdh*(BPG*Vmfor*nadph/(KBPG*Knadph) - Vmarev*gap*nadp*pi/(Kgap*Knadp*Kpi))/(Knadp*(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*pow(1 + nadph/Knadph + nadp/Knadp, 2)) - Vmarev*gap*pi*protGAPdh/(Kgap*Knadp*Kpi*(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*(1 + nadph/Knadph + nadp/Knadp));
    dwdx[11] = glc*kGDH;
    dwdx[12] = BPG*Vmfor*protGAPdh/(KBPG*Knadph*(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*(1 + nadph/Knadph + nadp/Knadp)) - protGAPdh*(BPG*Vmfor*nadph/(KBPG*Knadph) - Vmarev*gap*nadp*pi/(Kgap*Knadp*Kpi))/(Knadph*(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*pow(1 + nadph/Knadph + nadp/Knadp, 2));
    dwdx[13] = ADP*kPK;
    dwdx[14] = -protGAPdh*(1 + gap/Kgap)*(BPG*Vmfor*nadph/(KBPG*Knadph) - Vmarev*gap*nadp*pi/(Kgap*Knadp*Kpi))/(Kpi*pow(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi), 2)*(1 + nadph/Knadph + nadp/Knadp)) - Vmarev*gap*nadp*protGAPdh/(Kgap*Knadp*Kpi*(BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*(1 + nadph/Knadph + nadp/Knadp));
}