#include "sundials/sundials_types.h"

void JSparse_colptrs_kolodkin6(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 6;
    colptrs[3] = 10;
    colptrs[4] = 13;
    colptrs[5] = 19;
    colptrs[6] = 22;
    colptrs[7] = 25;
    colptrs[8] = 28;
}