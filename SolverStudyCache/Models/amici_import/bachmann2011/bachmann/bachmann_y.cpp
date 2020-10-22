#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_bachmann(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = CIS;
    y[1] = CISRNA;
    y[2] = CISnRNA1;
    y[3] = CISnRNA2;
    y[4] = CISnRNA3;
    y[5] = CISnRNA4;
    y[6] = CISnRNA5;
    y[7] = Epo;
    y[8] = EpoRJAK2;
    y[9] = EpoRJAK2_CIS;
    y[10] = EpoRpJAK2;
    y[11] = SHP1;
    y[12] = SHP1Act;
    y[13] = SOCS3;
    y[14] = SOCS3RNA;
    y[15] = SOCS3nRNA1;
    y[16] = SOCS3nRNA2;
    y[17] = SOCS3nRNA3;
    y[18] = SOCS3nRNA4;
    y[19] = SOCS3nRNA5;
    y[20] = STAT5;
    y[21] = npSTAT5;
    y[22] = p12EpoRpJAK2;
    y[23] = p1EpoRpJAK2;
    y[24] = p2EpoRpJAK2;
    y[25] = pSTAT5;
}