#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model1_kouril4(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 2;
    rowvals[3] = 3;
    rowvals[4] = 0;
    rowvals[5] = 1;
    rowvals[6] = 4;
    rowvals[7] = 5;
}