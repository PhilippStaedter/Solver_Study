#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_montanez1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*ARGin*Arginase_Vmaxarg/(ARGin + Arginase_Kmarg*(1 + ORN/Arginase_Kioarg));
    w[1] = 1.0*ARGex*Arginine_transport_Vmaxhat/((ARGex + Arginine_transport_Kmhat)*(ARGin/Arginine_transport_Kmhat + 1 + ORN/Arginine_transport_Kiornhat)) + 1.0*ARGex*Arginine_transport_Vmaxlat/((ARGex + Arginine_transport_Kmlat)*(ARGin/Arginine_transport_Kmlat + 1 + ORN/Arginine_transport_Kiornhat));
    w[2] = 1.0*ARGin*NOS_Vmaxnos1/(ARGin + NOS_Kmnos1);
    w[3] = 1.0*ODC_Vmaxodc*ORN/(ODC_Kmodc + ORN);
    w[4] = 1.0*ORN*Ornithine_efflux_Vmaxefflhat/((ORN + Ornithine_efflux_Kiornhat*(ARGin/Ornithine_efflux_Kmhat + 1))*(ARGex/Ornithine_efflux_Kmhat + 1)) + 1.0*ORN*Ornithine_efflux_Vmaxeffllat/((ORN + Ornithine_efflux_Kmeffllat*(ARGin/Ornithine_efflux_Kmlat + 1))*(ARGex/Ornithine_efflux_Kmlat + 1));
}