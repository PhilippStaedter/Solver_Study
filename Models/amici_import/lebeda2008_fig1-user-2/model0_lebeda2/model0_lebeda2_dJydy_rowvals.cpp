#include "sundials/sundials_types.h"

void dJydy_rowvals_model0_lebeda2(sunindextype *rowvals, int index){
    switch(index) {
        case 0:
                rowvals[0] = 0;
            break;
        case 1:
                rowvals[0] = 0;
            break;
        case 2:
                rowvals[0] = 0;
            break;
        case 3:
                rowvals[0] = 0;
            break;
    }
}