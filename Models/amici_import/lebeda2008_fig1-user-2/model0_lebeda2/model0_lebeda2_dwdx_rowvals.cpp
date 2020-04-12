#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_lebeda2(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 0;
    rowvals[2] = 2;
}