#include "sundials/sundials_types.h"

void dxdotdw_rowvals_kolodkin6(sunindextype *rowvals){
    rowvals[0] = 4;
    rowvals[1] = 6;
    rowvals[2] = 7;
    rowvals[3] = 1;
    rowvals[4] = 1;
    rowvals[5] = 2;
    rowvals[6] = 4;
    rowvals[7] = 3;
    rowvals[8] = 5;
    rowvals[9] = 2;
    rowvals[10] = 5;
    rowvals[11] = 3;
    rowvals[12] = 4;
}