#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_kouril5(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
}