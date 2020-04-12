#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"
#include "x.h"

void sx0_model3_levering2(realtype *sx0, const realtype t,const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 0:
            sx0[24] = 1;
            break;
        case 1:
            sx0[1] = 1;
            sx0[24] = -1;
            break;
    }
}