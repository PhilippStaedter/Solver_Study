#include "sundials/sundials_types.h"

void dxdotdw_rowvals_kolodkin2(sunindextype *rowvals){
    rowvals[0] = 2;
    rowvals[1] = 4;
    rowvals[2] = 5;
    rowvals[3] = 1;
    rowvals[4] = 1;
    rowvals[5] = 2;
    rowvals[6] = 3;
}