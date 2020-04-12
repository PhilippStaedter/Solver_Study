#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model1_kouril4(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -BPG*KeqPGK*VmPGK/(KmPGKATP*KmPGKP3G*(ADP/KmPGKADP + ATP/KmPGKATP + 1)*(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G)) + VmPGK*(ADP*BPG*KeqPGK - ATP*P3G)/(KmPGKADP*KmPGKATP*KmPGKP3G*pow(ADP/KmPGKADP + ATP/KmPGKATP + 1, 2)*(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G));
    dwdx[1] = kPK*pep;
    dwdx[2] = P3G*VmPGK/(KmPGKATP*KmPGKP3G*(ADP/KmPGKADP + ATP/KmPGKATP + 1)*(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G)) + VmPGK*(ADP*BPG*KeqPGK - ATP*P3G)/(pow(KmPGKATP, 2)*KmPGKP3G*pow(ADP/KmPGKADP + ATP/KmPGKATP + 1, 2)*(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G));
    dwdx[3] = -ADP*KeqPGK*VmPGK/(KmPGKATP*KmPGKP3G*(ADP/KmPGKADP + ATP/KmPGKATP + 1)*(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G)) + VmPGK*(ADP*BPG*KeqPGK - ATP*P3G)/(KmPGKATP*KmPGKBPG*KmPGKP3G*(ADP/KmPGKADP + ATP/KmPGKATP + 1)*pow(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G, 2));
    dwdx[4] = ATP*VmPGK/(KmPGKATP*KmPGKP3G*(ADP/KmPGKADP + ATP/KmPGKATP + 1)*(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G)) + VmPGK*(ADP*BPG*KeqPGK - ATP*P3G)/(KmPGKATP*pow(KmPGKP3G, 2)*(ADP/KmPGKADP + ATP/KmPGKATP + 1)*pow(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G, 2));
    dwdx[5] = ADP*kPK;
}