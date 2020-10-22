#include "sundials/sundials_types.h"

void JSparse_colptrs_kolodkin3(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 8;
    colptrs[3] = 11;
    colptrs[4] = 15;
    colptrs[5] = 20;
    colptrs[6] = 23;
    colptrs[7] = 27;
    colptrs[8] = 30;
    colptrs[9] = 33;
    colptrs[10] = 36;
}