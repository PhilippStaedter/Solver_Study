#include "sundials/sundials_types.h"

void JSparse_colptrs_kolodkin2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 4;
    colptrs[3] = 9;
    colptrs[4] = 12;
    colptrs[5] = 15;
    colptrs[6] = 18;
}