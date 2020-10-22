#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model2_band2(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 0;
    rowvals[2] = 0;
}