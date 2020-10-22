#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_bindschadler1(sunindextype *rowvals){
    rowvals[0] = 2;
    rowvals[1] = 3;
    rowvals[2] = 0;
    rowvals[3] = 1;
    rowvals[4] = 0;
    rowvals[5] = 1;
    rowvals[6] = 0;
    rowvals[7] = 1;
    rowvals[8] = 2;
    rowvals[9] = 3;
    rowvals[10] = 0;
    rowvals[11] = 1;
}