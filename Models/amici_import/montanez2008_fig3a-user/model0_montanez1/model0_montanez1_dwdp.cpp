#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_montanez1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*ARGin*Arginase_Kmarg*Arginase_Vmaxarg*ORN/(pow(Arginase_Kioarg, 2)*pow(ARGin + Arginase_Kmarg*(1 + ORN/Arginase_Kioarg), 2));
            break;
        case 1:
            dwdp[0] = 1.0*ARGin*Arginase_Vmaxarg*(-1 - ORN/Arginase_Kioarg)/pow(ARGin + Arginase_Kmarg*(1 + ORN/Arginase_Kioarg), 2);
            break;
        case 2:
            dwdp[0] = 1.0*ARGin/(ARGin + Arginase_Kmarg*(1 + ORN/Arginase_Kioarg));
            break;
        case 3:
            dwdp[1] = 1.0*ARGex*Arginine_transport_Vmaxhat*ORN/(pow(Arginine_transport_Kiornhat, 2)*(ARGex + Arginine_transport_Kmhat)*pow(ARGin/Arginine_transport_Kmhat + 1 + ORN/Arginine_transport_Kiornhat, 2)) + 1.0*ARGex*Arginine_transport_Vmaxlat*ORN/(pow(Arginine_transport_Kiornhat, 2)*(ARGex + Arginine_transport_Kmlat)*pow(ARGin/Arginine_transport_Kmlat + 1 + ORN/Arginine_transport_Kiornhat, 2));
            break;
        case 4:
            dwdp[1] = 1.0*ARGex/((ARGex + Arginine_transport_Kmlat)*(ARGin/Arginine_transport_Kmlat + 1 + ORN/Arginine_transport_Kiornhat));
            break;
        case 5:
            dwdp[1] = 1.0*ARGex*ARGin*Arginine_transport_Vmaxlat/(pow(Arginine_transport_Kmlat, 2)*(ARGex + Arginine_transport_Kmlat)*pow(ARGin/Arginine_transport_Kmlat + 1 + ORN/Arginine_transport_Kiornhat, 2)) - 1.0*ARGex*Arginine_transport_Vmaxlat/(pow(ARGex + Arginine_transport_Kmlat, 2)*(ARGin/Arginine_transport_Kmlat + 1 + ORN/Arginine_transport_Kiornhat));
            break;
        case 6:
            dwdp[1] = 1.0*ARGex/((ARGex + Arginine_transport_Kmhat)*(ARGin/Arginine_transport_Kmhat + 1 + ORN/Arginine_transport_Kiornhat));
            break;
        case 7:
            dwdp[1] = 1.0*ARGex*ARGin*Arginine_transport_Vmaxhat/(pow(Arginine_transport_Kmhat, 2)*(ARGex + Arginine_transport_Kmhat)*pow(ARGin/Arginine_transport_Kmhat + 1 + ORN/Arginine_transport_Kiornhat, 2)) - 1.0*ARGex*Arginine_transport_Vmaxhat/(pow(ARGex + Arginine_transport_Kmhat, 2)*(ARGin/Arginine_transport_Kmhat + 1 + ORN/Arginine_transport_Kiornhat));
            break;
        case 8:
            dwdp[2] = -1.0*ARGin*NOS_Vmaxnos1/pow(ARGin + NOS_Kmnos1, 2);
            break;
        case 9:
            dwdp[2] = 1.0*ARGin/(ARGin + NOS_Kmnos1);
            break;
        case 10:
            dwdp[3] = -1.0*ODC_Vmaxodc*ORN/pow(ODC_Kmodc + ORN, 2);
            break;
        case 11:
            dwdp[3] = 1.0*ORN/(ODC_Kmodc + ORN);
            break;
        case 12:
            dwdp[4] = 1.0*ORN*Ornithine_efflux_Vmaxefflhat*(-ARGin/Ornithine_efflux_Kmhat - 1)/(pow(ORN + Ornithine_efflux_Kiornhat*(ARGin/Ornithine_efflux_Kmhat + 1), 2)*(ARGex/Ornithine_efflux_Kmhat + 1));
            break;
        case 13:
            dwdp[4] = 1.0*ORN*Ornithine_efflux_Vmaxeffllat*(-ARGin/Ornithine_efflux_Kmlat - 1)/(pow(ORN + Ornithine_efflux_Kmeffllat*(ARGin/Ornithine_efflux_Kmlat + 1), 2)*(ARGex/Ornithine_efflux_Kmlat + 1));
            break;
        case 14:
            dwdp[4] = 1.0*ORN/((ORN + Ornithine_efflux_Kmeffllat*(ARGin/Ornithine_efflux_Kmlat + 1))*(ARGex/Ornithine_efflux_Kmlat + 1));
            break;
        case 15:
            dwdp[4] = 1.0*ORN/((ORN + Ornithine_efflux_Kiornhat*(ARGin/Ornithine_efflux_Kmhat + 1))*(ARGex/Ornithine_efflux_Kmhat + 1));
            break;
        case 16:
            dwdp[4] = 1.0*ARGex*ORN*Ornithine_efflux_Vmaxeffllat/(pow(Ornithine_efflux_Kmlat, 2)*(ORN + Ornithine_efflux_Kmeffllat*(ARGin/Ornithine_efflux_Kmlat + 1))*pow(ARGex/Ornithine_efflux_Kmlat + 1, 2)) + 1.0*ARGin*ORN*Ornithine_efflux_Kmeffllat*Ornithine_efflux_Vmaxeffllat/(pow(Ornithine_efflux_Kmlat, 2)*pow(ORN + Ornithine_efflux_Kmeffllat*(ARGin/Ornithine_efflux_Kmlat + 1), 2)*(ARGex/Ornithine_efflux_Kmlat + 1));
            break;
        case 17:
            dwdp[4] = 1.0*ARGex*ORN*Ornithine_efflux_Vmaxefflhat/(pow(Ornithine_efflux_Kmhat, 2)*(ORN + Ornithine_efflux_Kiornhat*(ARGin/Ornithine_efflux_Kmhat + 1))*pow(ARGex/Ornithine_efflux_Kmhat + 1, 2)) + 1.0*ARGin*ORN*Ornithine_efflux_Kiornhat*Ornithine_efflux_Vmaxefflhat/(pow(Ornithine_efflux_Kmhat, 2)*pow(ORN + Ornithine_efflux_Kiornhat*(ARGin/Ornithine_efflux_Kmhat + 1), 2)*(ARGex/Ornithine_efflux_Kmhat + 1));
            break;
    }
}