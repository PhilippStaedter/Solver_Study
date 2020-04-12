#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model24_karin3(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 2;
}