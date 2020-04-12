#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model22_karin1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 2;
    colptrs[3] = 3;
}