#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_montanez1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -1.0*ARGex*Arginine_transport_Vmaxhat/(pow(ARGex + Arginine_transport_Kmhat, 2)*(ARGin/Arginine_transport_Kmhat + 1 + ORN/Arginine_transport_Kiornhat)) - 1.0*ARGex*Arginine_transport_Vmaxlat/(pow(ARGex + Arginine_transport_Kmlat, 2)*(ARGin/Arginine_transport_Kmlat + 1 + ORN/Arginine_transport_Kiornhat)) + 1.0*Arginine_transport_Vmaxhat/((ARGex + Arginine_transport_Kmhat)*(ARGin/Arginine_transport_Kmhat + 1 + ORN/Arginine_transport_Kiornhat)) + 1.0*Arginine_transport_Vmaxlat/((ARGex + Arginine_transport_Kmlat)*(ARGin/Arginine_transport_Kmlat + 1 + ORN/Arginine_transport_Kiornhat));
    dwdx[1] = -1.0*ORN*Ornithine_efflux_Vmaxeffllat/(Ornithine_efflux_Kmlat*(ORN + Ornithine_efflux_Kmeffllat*(ARGin/Ornithine_efflux_Kmlat + 1))*pow(ARGex/Ornithine_efflux_Kmlat + 1, 2)) - 1.0*ORN*Ornithine_efflux_Vmaxefflhat/(Ornithine_efflux_Kmhat*(ORN + Ornithine_efflux_Kiornhat*(ARGin/Ornithine_efflux_Kmhat + 1))*pow(ARGex/Ornithine_efflux_Kmhat + 1, 2));
    dwdx[2] = -1.0*ARGin*Arginase_Vmaxarg/pow(ARGin + Arginase_Kmarg*(1 + ORN/Arginase_Kioarg), 2) + 1.0*Arginase_Vmaxarg/(ARGin + Arginase_Kmarg*(1 + ORN/Arginase_Kioarg));
    dwdx[3] = -1.0*ARGex*Arginine_transport_Vmaxlat/(Arginine_transport_Kmlat*(ARGex + Arginine_transport_Kmlat)*pow(ARGin/Arginine_transport_Kmlat + 1 + ORN/Arginine_transport_Kiornhat, 2)) - 1.0*ARGex*Arginine_transport_Vmaxhat/(Arginine_transport_Kmhat*(ARGex + Arginine_transport_Kmhat)*pow(ARGin/Arginine_transport_Kmhat + 1 + ORN/Arginine_transport_Kiornhat, 2));
    dwdx[4] = -1.0*ARGin*NOS_Vmaxnos1/pow(ARGin + NOS_Kmnos1, 2) + 1.0*NOS_Vmaxnos1/(ARGin + NOS_Kmnos1);
    dwdx[5] = -1.0*ORN*Ornithine_efflux_Kiornhat*Ornithine_efflux_Vmaxefflhat/(Ornithine_efflux_Kmhat*pow(ORN + Ornithine_efflux_Kiornhat*(ARGin/Ornithine_efflux_Kmhat + 1), 2)*(ARGex/Ornithine_efflux_Kmhat + 1)) - 1.0*ORN*Ornithine_efflux_Kmeffllat*Ornithine_efflux_Vmaxeffllat/(Ornithine_efflux_Kmlat*pow(ORN + Ornithine_efflux_Kmeffllat*(ARGin/Ornithine_efflux_Kmlat + 1), 2)*(ARGex/Ornithine_efflux_Kmlat + 1));
    dwdx[6] = -1.0*ARGin*Arginase_Kmarg*Arginase_Vmaxarg/(Arginase_Kioarg*pow(ARGin + Arginase_Kmarg*(1 + ORN/Arginase_Kioarg), 2));
    dwdx[7] = -1.0*ARGex*Arginine_transport_Vmaxhat/(Arginine_transport_Kiornhat*(ARGex + Arginine_transport_Kmhat)*pow(ARGin/Arginine_transport_Kmhat + 1 + ORN/Arginine_transport_Kiornhat, 2)) - 1.0*ARGex*Arginine_transport_Vmaxlat/(Arginine_transport_Kiornhat*(ARGex + Arginine_transport_Kmlat)*pow(ARGin/Arginine_transport_Kmlat + 1 + ORN/Arginine_transport_Kiornhat, 2));
    dwdx[8] = -1.0*ODC_Vmaxodc*ORN/pow(ODC_Kmodc + ORN, 2) + 1.0*ODC_Vmaxodc/(ODC_Kmodc + ORN);
    dwdx[9] = -1.0*ORN*Ornithine_efflux_Vmaxefflhat/(pow(ORN + Ornithine_efflux_Kiornhat*(ARGin/Ornithine_efflux_Kmhat + 1), 2)*(ARGex/Ornithine_efflux_Kmhat + 1)) - 1.0*ORN*Ornithine_efflux_Vmaxeffllat/(pow(ORN + Ornithine_efflux_Kmeffllat*(ARGin/Ornithine_efflux_Kmlat + 1), 2)*(ARGex/Ornithine_efflux_Kmlat + 1)) + 1.0*Ornithine_efflux_Vmaxefflhat/((ORN + Ornithine_efflux_Kiornhat*(ARGin/Ornithine_efflux_Kmhat + 1))*(ARGex/Ornithine_efflux_Kmhat + 1)) + 1.0*Ornithine_efflux_Vmaxeffllat/((ORN + Ornithine_efflux_Kmeffllat*(ARGin/Ornithine_efflux_Kmlat + 1))*(ARGex/Ornithine_efflux_Kmlat + 1));
}