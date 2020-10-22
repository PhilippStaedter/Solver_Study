#include "sundials/sundials_types.h"

void dJydy_colptrs_model2_(sunindextype *colptrs, int index){
    switch(index) {
        case 0:
                colptrs[0] = 0;
                colptrs[1] = 1;
                colptrs[2] = 1;
                colptrs[3] = 1;
                colptrs[4] = 1;
                colptrs[5] = 1;
                colptrs[6] = 1;
                colptrs[7] = 1;
            break;
        case 1:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 1;
                colptrs[3] = 1;
                colptrs[4] = 1;
                colptrs[5] = 1;
                colptrs[6] = 1;
                colptrs[7] = 1;
            break;
        case 2:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 0;
                colptrs[3] = 1;
                colptrs[4] = 1;
                colptrs[5] = 1;
                colptrs[6] = 1;
                colptrs[7] = 1;
            break;
        case 3:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 0;
                colptrs[3] = 0;
                colptrs[4] = 1;
                colptrs[5] = 1;
                colptrs[6] = 1;
                colptrs[7] = 1;
            break;
        case 4:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 0;
                colptrs[3] = 0;
                colptrs[4] = 0;
                colptrs[5] = 1;
                colptrs[6] = 1;
                colptrs[7] = 1;
            break;
        case 5:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 0;
                colptrs[3] = 0;
                colptrs[4] = 0;
                colptrs[5] = 0;
                colptrs[6] = 1;
                colptrs[7] = 1;
            break;
        case 6:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 0;
                colptrs[3] = 0;
                colptrs[4] = 0;
                colptrs[5] = 0;
                colptrs[6] = 0;
                colptrs[7] = 1;
            break;
    }
}