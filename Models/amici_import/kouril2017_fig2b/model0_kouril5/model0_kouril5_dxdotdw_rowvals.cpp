#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_kouril5(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
}