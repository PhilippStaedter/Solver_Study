#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_montanez1(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 2;
    rowvals[2] = 1;
    rowvals[3] = 1;
    rowvals[4] = 2;
    rowvals[5] = 2;
}