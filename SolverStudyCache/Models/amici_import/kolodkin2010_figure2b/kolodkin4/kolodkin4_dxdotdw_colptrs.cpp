#include "sundials/sundials_types.h"

void dxdotdw_colptrs_kolodkin4(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 4;
    colptrs[3] = 7;
    colptrs[4] = 9;
    colptrs[5] = 11;
    colptrs[6] = 13;
}