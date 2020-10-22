#include "sundials/sundials_types.h"

void dwdx_colptrs_kolodkin3(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 5;
    colptrs[3] = 7;
    colptrs[4] = 9;
    colptrs[5] = 11;
    colptrs[6] = 13;
    colptrs[7] = 15;
    colptrs[8] = 16;
    colptrs[9] = 17;
    colptrs[10] = 18;
}