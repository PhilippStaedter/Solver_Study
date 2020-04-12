#include "sundials/sundials_types.h"

void dxdotdw_rowvals_kolodkin3(sunindextype *rowvals){
    rowvals[0] = 4;
    rowvals[1] = 8;
    rowvals[2] = 9;
    rowvals[3] = 1;
    rowvals[4] = 1;
    rowvals[5] = 4;
    rowvals[6] = 7;
    rowvals[7] = 1;
    rowvals[8] = 3;
    rowvals[9] = 6;
    rowvals[10] = 5;
    rowvals[11] = 6;
    rowvals[12] = 2;
    rowvals[13] = 3;
    rowvals[14] = 2;
    rowvals[15] = 5;
}