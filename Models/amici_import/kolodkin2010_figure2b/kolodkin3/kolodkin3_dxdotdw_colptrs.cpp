#include "sundials/sundials_types.h"

void dxdotdw_colptrs_kolodkin3(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 4;
    colptrs[3] = 7;
    colptrs[4] = 10;
    colptrs[5] = 12;
    colptrs[6] = 14;
    colptrs[7] = 16;
}