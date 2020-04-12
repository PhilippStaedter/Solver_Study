#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model1_wang1(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 0;
    rowvals[2] = 0;
    rowvals[3] = 0;
    rowvals[4] = 0;
    rowvals[5] = 1;
    rowvals[6] = 1;
    rowvals[7] = 2;
    rowvals[8] = 2;
}