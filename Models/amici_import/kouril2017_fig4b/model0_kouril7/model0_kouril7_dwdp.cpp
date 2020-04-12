#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_kouril7(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -ADP*BPG*VmPGK/(KmPGKATP*KmPGKP3G*(ADP/KmPGKADP + ATP/KmPGKATP + 1)*(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G));
            break;
        case 1:
            dwdp[2] = BPG*(BPG*NADH*VmGAPDHr/(KmGAPDHBPG*KmGAPDHNADH) - NAD*VmGAPDHf*gap/(KmGAPDHGAP*KmGAPDHNAD))/(pow(KmGAPDHBPG, 2)*(1 + NADH/KmGAPDHNADH + NAD/KmGAPDHNAD)*pow(BPG/KmGAPDHBPG + 1 + gap/KmGAPDHGAP, 2)) - BPG*NADH*VmGAPDHr/(pow(KmGAPDHBPG, 2)*KmGAPDHNADH*(1 + NADH/KmGAPDHNADH + NAD/KmGAPDHNAD)*(BPG/KmGAPDHBPG + 1 + gap/KmGAPDHGAP));
            break;
        case 2:
            dwdp[2] = gap*(BPG*NADH*VmGAPDHr/(KmGAPDHBPG*KmGAPDHNADH) - NAD*VmGAPDHf*gap/(KmGAPDHGAP*KmGAPDHNAD))/(pow(KmGAPDHGAP, 2)*(1 + NADH/KmGAPDHNADH + NAD/KmGAPDHNAD)*pow(BPG/KmGAPDHBPG + 1 + gap/KmGAPDHGAP, 2)) + NAD*VmGAPDHf*gap/(pow(KmGAPDHGAP, 2)*KmGAPDHNAD*(1 + NADH/KmGAPDHNADH + NAD/KmGAPDHNAD)*(BPG/KmGAPDHBPG + 1 + gap/KmGAPDHGAP));
            break;
        case 3:
            dwdp[2] = NAD*(BPG*NADH*VmGAPDHr/(KmGAPDHBPG*KmGAPDHNADH) - NAD*VmGAPDHf*gap/(KmGAPDHGAP*KmGAPDHNAD))/(pow(KmGAPDHNAD, 2)*pow(1 + NADH/KmGAPDHNADH + NAD/KmGAPDHNAD, 2)*(BPG/KmGAPDHBPG + 1 + gap/KmGAPDHGAP)) + NAD*VmGAPDHf*gap/(KmGAPDHGAP*pow(KmGAPDHNAD, 2)*(1 + NADH/KmGAPDHNADH + NAD/KmGAPDHNAD)*(BPG/KmGAPDHBPG + 1 + gap/KmGAPDHGAP));
            break;
        case 4:
            dwdp[2] = -BPG*NADH*VmGAPDHr/(KmGAPDHBPG*pow(KmGAPDHNADH, 2)*(1 + NADH/KmGAPDHNADH + NAD/KmGAPDHNAD)*(BPG/KmGAPDHBPG + 1 + gap/KmGAPDHGAP)) + NADH*(BPG*NADH*VmGAPDHr/(KmGAPDHBPG*KmGAPDHNADH) - NAD*VmGAPDHf*gap/(KmGAPDHGAP*KmGAPDHNAD))/(pow(KmGAPDHNADH, 2)*pow(1 + NADH/KmGAPDHNADH + NAD/KmGAPDHNAD, 2)*(BPG/KmGAPDHBPG + 1 + gap/KmGAPDHGAP));
            break;
        case 5:
            dwdp[0] = -ADP*VmPGK*(ADP*BPG*KeqPGK - ATP*P3G)/(pow(KmPGKADP, 2)*KmPGKATP*KmPGKP3G*pow(ADP/KmPGKADP + ATP/KmPGKATP + 1, 2)*(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G));
            break;
        case 6:
            dwdp[0] = -ATP*VmPGK*(ADP*BPG*KeqPGK - ATP*P3G)/(pow(KmPGKATP, 3)*KmPGKP3G*pow(ADP/KmPGKADP + ATP/KmPGKATP + 1, 2)*(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G)) + VmPGK*(ADP*BPG*KeqPGK - ATP*P3G)/(pow(KmPGKATP, 2)*KmPGKP3G*(ADP/KmPGKADP + ATP/KmPGKATP + 1)*(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G));
            break;
        case 7:
            dwdp[0] = -BPG*VmPGK*(ADP*BPG*KeqPGK - ATP*P3G)/(KmPGKATP*pow(KmPGKBPG, 2)*KmPGKP3G*(ADP/KmPGKADP + ATP/KmPGKATP + 1)*pow(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G, 2));
            break;
        case 8:
            dwdp[0] = VmPGK*(ADP*BPG*KeqPGK - ATP*P3G)/(KmPGKATP*pow(KmPGKP3G, 2)*(ADP/KmPGKADP + ATP/KmPGKATP + 1)*(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G)) - P3G*VmPGK*(ADP*BPG*KeqPGK - ATP*P3G)/(KmPGKATP*pow(KmPGKP3G, 3)*(ADP/KmPGKADP + ATP/KmPGKATP + 1)*pow(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G, 2));
            break;
        case 9:
            dwdp[2] = -NAD*gap/(KmGAPDHGAP*KmGAPDHNAD*(1 + NADH/KmGAPDHNADH + NAD/KmGAPDHNAD)*(BPG/KmGAPDHBPG + 1 + gap/KmGAPDHGAP));
            break;
        case 10:
            dwdp[2] = BPG*NADH/(KmGAPDHBPG*KmGAPDHNADH*(1 + NADH/KmGAPDHNADH + NAD/KmGAPDHNAD)*(BPG/KmGAPDHBPG + 1 + gap/KmGAPDHGAP));
            break;
        case 11:
            dwdp[0] = -(ADP*BPG*KeqPGK - ATP*P3G)/(KmPGKATP*KmPGKP3G*(ADP/KmPGKADP + ATP/KmPGKATP + 1)*(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G));
            break;
        case 12:
            dwdp[3] = Glc*NAD;
            break;
        case 13:
            dwdp[1] = ADP*pep;
            break;
    }
}