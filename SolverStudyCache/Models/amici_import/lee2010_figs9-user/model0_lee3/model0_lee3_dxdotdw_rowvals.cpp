#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_lee3(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 2;
    rowvals[2] = 1;
    rowvals[3] = 2;
    rowvals[4] = 0;
    rowvals[5] = 3;
    rowvals[6] = 1;
    rowvals[7] = 3;
}