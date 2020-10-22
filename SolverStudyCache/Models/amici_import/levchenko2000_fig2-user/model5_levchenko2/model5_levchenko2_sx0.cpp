#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"
#include "x.h"

void sx0_model5_levchenko2(realtype *sx0, const realtype t,const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 0:
            sx0[0] = initial_value_Total_Scaffold;
            break;
        case 1:
            sx0[1] = initial_value_Total_Scaffold;
            break;
        case 2:
            sx0[3] = initial_value_Total_Scaffold;
            break;
        case 3:
            sx0[5] = initial_value_Total_Scaffold;
            break;
        case 4:
            sx0[9] = 1;
            break;
        case 5:
            sx0[17] = 1;
            break;
        case 6:
            sx0[26] = 1;
            break;
        case 7:
            sx0[0] = initial_value_C1_Fraction;
            sx0[1] = initial_value_C2_Fraction;
            sx0[3] = initial_value_C4_Fraction;
            sx0[5] = initial_value_C6_Fraction;
            break;
    }
}