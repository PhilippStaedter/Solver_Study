#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_kouril6(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -BPG*VmrPGK*protPGK/(KpgkADP*KpgkBPG*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G)) - BPG*protPGK*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(KpgkADP*KpgkBPG*(ADP/KiADP + 1)*pow(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G, 2)) - protPGK*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(KiADP*pow(ADP/KiADP + 1, 2)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G));
    dwdx[1] = kPK*pep;
    dwdx[2] = P3G*VmfPGK*protPGK/(KpgkATP*KpgkP3G*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G)) - P3G*protPGK*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(KpgkATP*KpgkP3G*(ADP/KiADP + 1)*pow(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G, 2));
    dwdx[3] = -ADP*VmrPGK*protPGK/(KpgkADP*KpgkBPG*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G)) - protPGK*(ADP/KpgkADP + 1)*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(KpgkBPG*(ADP/KiADP + 1)*pow(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G, 2));
    dwdx[4] = kd;
    dwdx[5] = ATP*VmfPGK*protPGK/(KpgkATP*KpgkP3G*(ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G)) - protPGK*(ATP/KpgkATP + 1)*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/(KpgkP3G*(ADP/KiADP + 1)*pow(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G, 2));
    dwdx[6] = ADP*kPK;
}