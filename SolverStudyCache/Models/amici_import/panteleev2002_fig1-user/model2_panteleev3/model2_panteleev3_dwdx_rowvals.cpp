#include "sundials/sundials_types.h"

void dwdx_rowvals_model2_panteleev3(sunindextype *rowvals){
    rowvals[0] = 3;
    rowvals[1] = 0;
    rowvals[2] = 2;
    rowvals[3] = 4;
    rowvals[4] = 0;
    rowvals[5] = 1;
    rowvals[6] = 2;
    rowvals[7] = 0;
    rowvals[8] = 2;
    rowvals[9] = 3;
    rowvals[10] = 3;
    rowvals[11] = 4;
    rowvals[12] = 4;
}