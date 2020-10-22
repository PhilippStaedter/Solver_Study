#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_kouril6(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 2;
    rowvals[3] = 3;
    rowvals[4] = 0;
    rowvals[5] = 1;
    rowvals[6] = 4;
    rowvals[7] = 6;
    rowvals[8] = 2;
    rowvals[9] = 3;
    rowvals[10] = 5;
}