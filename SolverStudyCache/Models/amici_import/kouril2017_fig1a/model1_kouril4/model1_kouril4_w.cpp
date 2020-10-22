#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model1_kouril4(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = -VmPGK*(ADP*BPG*KeqPGK - ATP*P3G)/(KmPGKATP*KmPGKP3G*(ADP/KmPGKADP + ATP/KmPGKATP + 1)*(BPG/KmPGKBPG + 1 + P3G/KmPGKP3G));
    w[1] = ADP*kPK*pep;
}