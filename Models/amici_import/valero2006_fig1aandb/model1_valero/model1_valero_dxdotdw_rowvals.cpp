#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model1_valero(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 2;
    rowvals[2] = 0;
    rowvals[3] = 1;
    rowvals[4] = 2;
    rowvals[5] = 0;
    rowvals[6] = 2;
    rowvals[7] = 5;
    rowvals[8] = 4;
    rowvals[9] = 5;
}