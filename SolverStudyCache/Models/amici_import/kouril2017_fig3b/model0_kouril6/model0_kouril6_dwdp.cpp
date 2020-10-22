#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_kouril6(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = ADP*protPGK*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(pow(KiADP, 2)*pow(ADP/KiADP + 1, 2)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G));
            break;
        case 1:
            dwdp[0] = ADP*BPG*VmrPGK*protPGK/(pow(KpgkADP, 2)*KpgkBPG*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G)) + ADP*BPG*protPGK*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(pow(KpgkADP, 2)*KpgkBPG*(ADP/KiADP + 1)*pow(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G, 2));
            break;
        case 2:
            dwdp[0] = -ATP*P3G*VmfPGK*protPGK/(pow(KpgkATP, 2)*KpgkP3G*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G)) + ATP*P3G*protPGK*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(pow(KpgkATP, 2)*KpgkP3G*(ADP/KiADP + 1)*pow(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G, 2));
            break;
        case 3:
            dwdp[0] = ADP*BPG*VmrPGK*protPGK/(KpgkADP*pow(KpgkBPG, 2)*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G)) + BPG*protPGK*(ADP/KpgkADP + 1)*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(pow(KpgkBPG, 2)*(ADP/KiADP + 1)*pow(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G, 2));
            break;
        case 4:
            dwdp[0] = -ATP*P3G*VmfPGK*protPGK/(KpgkATP*pow(KpgkP3G, 2)*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G)) + P3G*protPGK*(ATP/KpgkATP + 1)*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(pow(KpgkP3G, 2)*(ADP/KiADP + 1)*pow(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G, 2));
            break;
        case 5:
            dwdp[0] = ATP*P3G*protPGK/(KpgkATP*KpgkP3G*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G));
            break;
        case 6:
            dwdp[0] = -ADP*BPG*protPGK/(KpgkADP*KpgkBPG*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G));
            break;
        case 7:
            dwdp[1] = ADP*pep;
            break;
        case 8:
            dwdp[2] = BPG;
            break;
        case 9:
            dwdp[0] = (-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/((ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G));
            break;
    }
}