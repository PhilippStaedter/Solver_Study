#include "sundials/sundials_types.h"

void dwdx_rowvals_model2_band2(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
}