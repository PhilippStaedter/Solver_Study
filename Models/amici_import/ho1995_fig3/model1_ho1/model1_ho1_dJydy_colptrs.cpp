#include "sundials/sundials_types.h"

void dJydy_colptrs_model1_ho1(sunindextype *colptrs, int index){
    switch(index) {
        case 0:
                colptrs[0] = 0;
                colptrs[1] = 1;
            break;
    }
}