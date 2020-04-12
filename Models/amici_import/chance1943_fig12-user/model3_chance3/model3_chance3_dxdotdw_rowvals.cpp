#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model3_chance3(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 3;
    rowvals[3] = 0;
    rowvals[4] = 1;
    rowvals[5] = 2;
}