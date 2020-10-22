#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model6_panteleev3(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 2;
    rowvals[2] = 4;
    rowvals[3] = 2;
    rowvals[4] = 3;
    rowvals[5] = 1;
    rowvals[6] = 3;
    rowvals[7] = 5;
    rowvals[8] = 0;
    rowvals[9] = 5;
    rowvals[10] = 6;
    rowvals[11] = 1;
    rowvals[12] = 6;
    rowvals[13] = 7;
}