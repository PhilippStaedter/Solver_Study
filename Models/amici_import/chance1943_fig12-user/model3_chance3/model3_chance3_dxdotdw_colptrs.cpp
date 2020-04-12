#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model3_chance3(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 6;
}