#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_kouril8(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = protPGK*(-ADP*BPG*VmrPGK/(KpgkADP*KpgkBPG) + ATP*P3G*VmfPGK/(KpgkATP*KpgkP3G))/((ADP/KiADP + 1)*(BPG*(ADP/KpgkADP + 1)/KpgkBPG + 1 + P3G*(ATP/KpgkATP + 1)/KpgkP3G));
    w[1] = ADP*kPK*pep;
    w[2] = BPG*kdbpg;
    w[3] = protGAPdh*(BPG*Vmfor*nadph/(KBPG*Knadph) - Vmarev*gap*nadp*pi/(Kgap*Knadp*Kpi))/((BPG/KBPG + (1 + gap/Kgap)*(1 + pi/Kpi))*(1 + nadph/Knadph + nadp/Knadp));
    w[4] = glc*kGDH*nadp;
    w[5] = gap*kdgap;
}