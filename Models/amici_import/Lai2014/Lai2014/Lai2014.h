#ifndef _amici_TPL_MODELNAME_h
#define _amici_TPL_MODELNAME_h
#include <cmath>
#include <memory>

#include "amici/model_ode.h"
#include "amici/solver_cvodes.h"

#include "sundials/sundials_types.h"

namespace amici {
class Solver;
}

/**
 * @brief Wrapper function to instantiate the linked Amici model without knowing
 * the name at compile time.
 * @return
 */
extern void J_Lai2014(realtype *J, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w,
                            const realtype *dwdx);
extern void JB_Lai2014(realtype *JB, const realtype t, const realtype *x,
                             const realtype *p, const realtype *k,
                             const realtype *h, const realtype *xB,
                             const realtype *w, const realtype *dwdx);
extern void JDiag_Lai2014(realtype *JDiag, const realtype t,
                                const realtype *x, const realtype *p,
                                const realtype *k, const realtype *h,
                                const realtype *w, const realtype *dwdx);
extern void JSparse_Lai2014(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JSparse_colptrs_Lai2014(sunindextype *colptrs);
extern void JSparse_rowvals_Lai2014(sunindextype *rowvals);
extern void JSparseB_Lai2014(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void JSparseB_colptrs_Lai2014(sunindextype *colptrs);
extern void JSparseB_rowvals_Lai2014(sunindextype *rowvals);
extern void Jy_Lai2014(realtype *nllh, const int iy, const realtype *p,
                             const realtype *k, const realtype *y,
                             const realtype *sigmay, const realtype *my);
extern void dJydsigmay_Lai2014(realtype *dJydsigmay, const int iy,
                                     const realtype *p, const realtype *k,
                                     const realtype *y, const realtype *sigmay,
                                     const realtype *my);
extern void dJydy_Lai2014(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_Lai2014(sunindextype *colptrs, int index);
extern void dJydy_rowvals_Lai2014(sunindextype *rowvals, int index);
extern void dwdp_Lai2014(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip);
extern void dwdx_Lai2014(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl);
extern void dwdx_colptrs_Lai2014(sunindextype *colptrs);
extern void dwdx_rowvals_Lai2014(sunindextype *rowvals);
extern void dxdotdw_Lai2014(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdw_colptrs_Lai2014(sunindextype *colptrs);
extern void dxdotdw_rowvals_Lai2014(sunindextype *rowvals);
extern void dxdotdp_Lai2014(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w);
extern void dydx_Lai2014(realtype *dydx, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w, const realtype *dwdx);
extern void dydp_Lai2014(realtype *dydp, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const int ip, const realtype *w,
                               const realtype *dwp);
extern void dsigmaydp_Lai2014(realtype *dsigmaydp, const realtype t,
                                    const realtype *p, const realtype *k,
                                    const int ip);
extern void sigmay_Lai2014(realtype *sigmay, const realtype t,
                                 const realtype *p, const realtype *k);
extern void w_Lai2014(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl);
extern void x0_Lai2014(realtype *x0, const realtype t, const realtype *p,
                             const realtype *k);
extern void x0_fixedParameters_Lai2014(realtype *x0, const realtype t,
                                             const realtype *p,
                                             const realtype *k);
extern void sx0_Lai2014(realtype *sx0, const realtype t,
                              const realtype *x0, const realtype *p,
                              const realtype *k, const int ip);
extern void sx0_fixedParameters_Lai2014(realtype *sx0, const realtype t,
                                              const realtype *x0,
                                              const realtype *p,
                                              const realtype *k, const int ip);
extern void xdot_Lai2014(realtype *xdot, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w);
extern void y_Lai2014(realtype *y, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w);

extern void x_solver_Lai2014(realtype *x_solver, const realtype *x_rdata);
extern void total_cl_Lai2014(realtype *total_cl, const realtype *x_rdata);

/**
 * @brief AMICI-generated model subclass.
 */
class Model_Lai2014 : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_Lai2014()
        : amici::Model_ODE(
              195,                                // nx_rdata
              195,                            // nxtrue_rdata
              195,                               // nx_solver
              195,                           // nxtrue_solver
              195,                                      // ny
              195,                                  // nytrue
              0,                                      // nz
              0,                                  // nztrue
              0,                                  // nevent
              1,                              // nobjective
              576,                                      // nw
              1664,                                   // ndwdx
              15552,                                   // ndwdp
              1280,                                // ndxdotdw
              std::vector<int>{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},                                  // ndjydy
              2050,                                     // nnz
              195,                                     // ubw
              195,                                     // lbw
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{8616.61, 322686.0, 0.0006, 100000000.0, 0.0006, 1000000000.0, 1000000000.0, 10000000.0, 0.000215, 10000.0, 0.000317, 10000000.0, 10000.0, 1.0, 6.242e-05, 1e-06, 9.192e-05, 5e-11, 7e-08, 1.0, 1.0, 1e-06, 1e-09, 1e-06, 100000000.0, 6.242e-05, 9.192e-05},       // dynamic parameters
              std::vector<realtype>{}, // fixedParameters
              std::vector<int>{},                          // plist
              std::vector<realtype>(195, 0.0),   // idlist
              std::vector<int>{}                           // z2event
          ) {}

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    virtual amici::Model *clone() const override {
        return new Model_Lai2014(*this);
    }

    /** model specific implementation for fJ
     * @param J Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJ(realtype *J, const realtype t, const realtype *x,
                    const realtype *p, const realtype *k, const realtype *h,
                    const realtype *w, const realtype *dwdx) override {
        J_Lai2014(J, t, x, p, k, h, w, dwdx);
    }

    /** model specific implementation for fJB
     * @param JB Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param xB Vector with the adjoint states
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJB(realtype *JB, const realtype t, const realtype *x,
                     const realtype *p, const realtype *k, const realtype *h,
                     const realtype *xB, const realtype *w,
                     const realtype *dwdx) override {
        JB_Lai2014(JB, t, x, p, k, h, xB, w, dwdx);
    }

    /** model specific implementation for fJDiag
     * @param JDiag Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJDiag(realtype *JDiag, const realtype t, const realtype *x,
                        const realtype *p, const realtype *k, const realtype *h,
                        const realtype *w, const realtype *dwdx) override {
        JDiag_Lai2014(JDiag, t, x, p, k, h, w, dwdx);
    }

    virtual void fJSparse(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        JSparse_Lai2014( JSparse,  t,  x,  p,  k,  h,  w,  dwdx);
    }


    virtual void fJSparse_colptrs(sunindextype *colptrs) override {
        JSparse_colptrs_Lai2014(colptrs);
    }


    virtual void fJSparse_rowvals(sunindextype *rowvals) override {
        JSparse_rowvals_Lai2014(rowvals);
    }


    virtual void fJSparseB(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        JSparseB_Lai2014( JSparseB,  t,  x,  p,  k,  h,  xB,  w,  dwdx);
    }


    virtual void fJSparseB_colptrs(sunindextype *colptrs) override {
        JSparseB_colptrs_Lai2014(colptrs);
    }


    virtual void fJSparseB_rowvals(sunindextype *rowvals) override {
        JSparseB_rowvals_Lai2014(rowvals);
    }


    /** model specific implementation of fJrz
     * @param nllh regularization for event measurements z
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     **/
    virtual void fJrz(realtype *nllh, const int iz, const realtype *p,
                      const realtype *k, const realtype *rz,
                      const realtype *sigmaz) override {}

    /** model specific implementation of fJy
     * @param nllh negative log-likelihood for measurements y
     * @param iy output index
     * @param p parameter vector
     * @param k constant vector
     * @param y model output at timepoint
     * @param sigmay measurement standard deviation at timepoint
     * @param my measurements at timepoint
     **/
    virtual void fJy(realtype *nllh, const int iy, const realtype *p,
                     const realtype *k, const realtype *y,
                     const realtype *sigmay, const realtype *my) override {
        Jy_Lai2014(nllh, iy, p, k, y, sigmay, my);
    }

    /** model specific implementation of fJz
     * @param nllh negative log-likelihood for event measurements z
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurements at timepoint
     **/
    virtual void fJz(realtype *nllh, const int iz, const realtype *p,
                     const realtype *k, const realtype *z,
                     const realtype *sigmaz, const realtype *mz) override {}

    /** model specific implementation of fdJrzdsigma
     * @param dJrzdsigma Sensitivity of event penalization Jrz w.r.t.
     * standard deviation sigmaz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param rz model root output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     **/
    virtual void fdJrzdsigma(realtype *dJrzdsigma, const int iz,
                             const realtype *p, const realtype *k,
                             const realtype *rz,
                             const realtype *sigmaz) override {}

    /** model specific implementation of fdJrzdz
     * @param dJrzdz partial derivative of event penalization Jrz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param rz model root output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     **/
    virtual void fdJrzdz(realtype *dJrzdz, const int iz, const realtype *p,
                         const realtype *k, const realtype *rz,
                         const realtype *sigmaz) override {}

    /** model specific implementation of fdJydsigma
     * @param dJydsigma Sensitivity of time-resolved measurement
     * negative log-likelihood Jy w.r.t. standard deviation sigmay
     * @param iy output index
     * @param p parameter vector
     * @param k constant vector
     * @param y model output at timepoint
     * @param sigmay measurement standard deviation at timepoint
     * @param my measurement at timepoint
     **/
    virtual void fdJydsigma(realtype *dJydsigma, const int iy,
                            const realtype *p, const realtype *k,
                            const realtype *y, const realtype *sigmay,
                            const realtype *my) override {
        dJydsigmay_Lai2014(dJydsigma, iy, p, k, y, sigmay, my);
    }


    /** model specific implementation of fdJzdsigma
     * @param dJzdsigma Sensitivity of event measurement
     * negative log-likelihood Jz w.r.t. standard deviation sigmaz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurement at timepoint
     **/
    virtual void fdJzdsigma(realtype *dJzdsigma, const int iz,
                            const realtype *p, const realtype *k,
                            const realtype *z, const realtype *sigmaz,
                            const realtype *mz) override {}

    /** model specific implementation of fdJzdz
     * @param dJzdz partial derivative of event measurement negative
     *log-likelihood Jz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurement at timepoint
     **/
    virtual void fdJzdz(realtype *dJzdz, const int iz, const realtype *p,
                        const realtype *k, const realtype *z,
                        const realtype *sigmaz, const realtype *mz) override {}

    /** model specific implementation of fdeltasx
     * @param deltaqB sensitivity update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ip sensitivity index
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param xB adjoint state
     **/
    virtual void fdeltaqB(realtype *deltaqB, const realtype t,
                          const realtype *x, const realtype *p,
                          const realtype *k, const realtype *h, const int ip,
                          const int ie, const realtype *xdot,
                          const realtype *xdot_old,
                          const realtype *xB) override {}

    /** model specific implementation of fdeltasx
     * @param deltasx sensitivity update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param w repeating elements vector
     * @param ip sensitivity index
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param sx state sensitivity
     * @param stau event-time sensitivity
     **/
    virtual void fdeltasx(realtype *deltasx, const realtype t,
                          const realtype *x, const realtype *p,
                          const realtype *k, const realtype *h,
                          const realtype *w, const int ip, const int ie,
                          const realtype *xdot, const realtype *xdot_old,
                          const realtype *sx, const realtype *stau) override {}

    /** model specific implementation of fdeltax
     * @param deltax state update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     **/
    virtual void fdeltax(realtype *deltax, const realtype t, const realtype *x,
                         const realtype *p, const realtype *k,
                         const realtype *h, const int ie, const realtype *xdot,
                         const realtype *xdot_old) override {}

    /** model specific implementation of fdeltaxB
     * @param deltaxB adjoint state update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param xB current adjoint state
     **/
    virtual void fdeltaxB(realtype *deltaxB, const realtype t,
                          const realtype *x, const realtype *p,
                          const realtype *k, const realtype *h, const int ie,
                          const realtype *xdot, const realtype *xdot_old,
                          const realtype *xB) override {}

    /** model specific implementation of fdrzdp
     * @param drzdp partial derivative of root output rz w.r.t. model parameters
     *p
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     **/
    virtual void fdrzdp(realtype *drzdp, const int ie, const realtype t,
                        const realtype *x, const realtype *p, const realtype *k,
                        const realtype *h, const int ip) override {}

    /** model specific implementation of fdrzdx
     * @param drzdx partial derivative of root output rz w.r.t. model states x
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fdrzdx(realtype *drzdx, const int ie, const realtype t,
                        const realtype *x, const realtype *p, const realtype *k,
                        const realtype *h) override {}

    /** model specific implementation of fsigmay
     * @param dsigmaydp partial derivative of standard deviation of measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fdsigmaydp(realtype *dsigmaydp, const realtype t,
                            const realtype *p, const realtype *k,
                            const int ip) override {
        dsigmaydp_Lai2014(dsigmaydp, t, p, k, ip);
    }

    /** model specific implementation of fsigmaz
     * @param dsigmazdp partial derivative of standard deviation of event
     *measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fdsigmazdp(realtype *dsigmazdp, const realtype t,
                            const realtype *p, const realtype *k,
                            const int ip) override {}

    virtual void fdJydy(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydy_Lai2014( dJydy,  iy,  p,  k,  y,  sigmay,  my);
    }

    virtual void fdJydy_colptrs(sunindextype *colptrs, int index) override {
        dJydy_colptrs_Lai2014(colptrs, index);
    }

    virtual void fdJydy_rowvals(sunindextype *rowvals, int index) override {
        dJydy_rowvals_Lai2014(rowvals, index);
    }


    virtual void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip) override {
        dwdp_Lai2014( dwdp,  t,  x,  p,  k,  h,  w,  tcl,  dtcldp,  ip);
    }


    virtual void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {
        dwdx_Lai2014( dwdx,  t,  x,  p,  k,  h,  w,  tcl);
    }


    virtual void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdw_Lai2014( dxdotdw,  t,  x,  p,  k,  h,  w);
    }


    virtual void fdxdotdw_colptrs(sunindextype *colptrs) override {
        dxdotdw_colptrs_Lai2014(colptrs);
    }


    virtual void fdxdotdw_rowvals(sunindextype *rowvals) override {
        dxdotdw_rowvals_Lai2014(rowvals);
    }


    virtual void fdxdotdp(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w) override {
        dxdotdp_Lai2014( dxdotdp,  t,  x,  p,  k,  h,  ip,  w);
    }


    /** model specific implementation of fdydx
     * @param dydx partial derivative of observables y w.r.t. model states x
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fdydx(realtype *dydx, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const realtype *w, const realtype *dwdx) override {
        dydx_Lai2014(dydx, t, x, p, k, h, w, dwdx);
    }

    /** model specific implementation of fdydp
     * @param dydp partial derivative of observables y w.r.t. model parameters p
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     **/
    virtual void fdydp(realtype *dydp, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const int ip, const realtype *w,
                       const realtype *dwdp) override {
        dydp_Lai2014(dydp, t, x, p, k, h, ip, w, dwdp);
    }

    /** model specific implementation of fdzdp
     * @param dzdp partial derivative of event-resolved output z w.r.t. model
     *parameters p
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     **/
    virtual void fdzdp(realtype *dzdp, const int ie, const realtype t,
                       const realtype *x, const realtype *p, const realtype *k,
                       const realtype *h, const int ip) override {}

    /** model specific implementation of fdzdx
     * @param dzdx partial derivative of event-resolved output z w.r.t. model
     *states x
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fdzdx(realtype *dzdx, const int ie, const realtype t,
                       const realtype *x, const realtype *p, const realtype *k,
                       const realtype *h) override {}

    /** model specific implementation for froot
     * @param root values of the trigger function
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     **/
    virtual void froot(realtype *root, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k,
                       const realtype *h) override {}

    /** model specific implementation of frz
     * @param rz value of root function at current timepoint (non-output events
     *not included)
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void frz(realtype *rz, const int ie, const realtype t,
                     const realtype *x, const realtype *p, const realtype *k,
                     const realtype *h) override {}

    /** model specific implementation of fsigmay
     * @param sigmay standard deviation of measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fsigmay(realtype *sigmay, const realtype t, const realtype *p,
                         const realtype *k) override {
        sigmay_Lai2014(sigmay, t, p, k);
    }

    /** model specific implementation of fsigmaz
     * @param sigmaz standard deviation of event measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fsigmaz(realtype *sigmaz, const realtype t, const realtype *p,
                         const realtype *k) override {}

    /** model specific implementation of fsrz
     * @param srz Sensitivity of rz, total derivative
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param sx current state sensitivity
     * @param h heavyside vector
     * @param ip sensitivity index
     **/
    virtual void fsrz(realtype *srz, const int ie, const realtype t,
                      const realtype *x, const realtype *p, const realtype *k,
                      const realtype *h, const realtype *sx,
                      const int ip) override {}

    /** model specific implementation of fstau
     * @param stau total derivative of event timepoint
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param sx current state sensitivity
     * @param ip sensitivity index
     * @param ie event index
     **/
    virtual void fstau(realtype *stau, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const realtype *sx, const int ip,
                       const int ie) override {}

    /** model specific implementation of fsx0
     * @param sx0 initial state sensitivities
     * @param t initial time
     * @param x0 initial state
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fsx0(realtype *sx0, const realtype t, const realtype *x0,
                      const realtype *p, const realtype *k,
                      const int ip) override {
        sx0_Lai2014(sx0, t, x0, p, k, ip);
    }

    /** model specific implementation of fsx0_fixedParameters
     * @param sx0 initial state sensitivities
     * @param t initial time
     * @param x0 initial state
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fsx0_fixedParameters(realtype *sx0, const realtype t,
                                      const realtype *x0, const realtype *p,
                                      const realtype *k,
                                      const int ip) override {
        sx0_fixedParameters_Lai2014(sx0, t, x0, p, k, ip);
    }

    /** model specific implementation of fsz
     * @param sz Sensitivity of rz, total derivative
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param sx current state sensitivity
     * @param ip sensitivity index
     **/
    virtual void fsz(realtype *sz, const int ie, const realtype t,
                     const realtype *x, const realtype *p, const realtype *k,
                     const realtype *h, const realtype *sx,
                     const int ip) override {}

    virtual void fw(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl) override {
        w_Lai2014( w,  t,  x,  p,  k,  h,  tcl);
    }


    /** model specific implementation of fx0
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0(realtype *x0, const realtype t, const realtype *p,
                     const realtype *k) override {
        x0_Lai2014(x0, t, p, k);
    }

    /** model specific implementation of fx0_fixedParameters
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0_fixedParameters(realtype *x0, const realtype t,
                                     const realtype *p,
                                     const realtype *k) override {
        x0_fixedParameters_Lai2014(x0, t, p, k);
    }

    /** model specific implementation for fxdot
     * @param xdot residual function
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     **/
    virtual void fxdot(realtype *xdot, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const realtype *w) override {
        xdot_Lai2014(xdot, t, x, p, k, h, w);
    }

    /** model specific implementation of fy
     * @param y model output at current timepoint
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fy(realtype *y, const realtype t, const realtype *x,
                    const realtype *p, const realtype *k, const realtype *h,
                    const realtype *w) override {
        y_Lai2014(y, t, x, p, k, h, w);
    }

    /** model specific implementation of fz
     * @param z value of event output
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fz(realtype *z, const int ie, const realtype t,
                    const realtype *x, const realtype *p, const realtype *k,
                    const realtype *h) override {}

    

    virtual void fx_solver(realtype *x_solver, const realtype *x_rdata) override {
        x_solver_Lai2014( x_solver,  x_rdata);
    }


    virtual void ftotal_cl(realtype *total_cl, const realtype *x_rdata) override {
        total_cl_Lai2014( total_cl,  x_rdata);
    }


    /**
     * @brief Get names of the model parameters
     * @return the names
     */
    virtual std::vector<std::string> getParameterNames() const override {
        return std::vector<std::string>{"lC",
"lN",
"Kd_rbp_TT",
"kon_tbp",
"Kd_rbp_RT",
"kon_AT",
"kon_AR",
"kon_CR",
"cN",
"k_R2T_C",
"cC",
"kon_CT",
"k_R2T_N",
"Kd_tbp_RT",
"KDT",
"conc_rbp",
"KBT",
"Kd_rbp_RR",
"Kd_rbp_TR",
"Kd_tbp_TT",
"Kd_tbp_TR",
"conc_cam",
"Kd_tbp_RR",
"conc_tbp",
"kon_rbp",
"KCT",
"KAT",};
    }

    /**
     * @brief Get names of the model states
     * @return the names
     */
    virtual std::vector<std::string> getStateNames() const override {
        return std::vector<std::string>{"cam_RR_0_0",
"cam_RR_0_rbp",
"cam_RR_0_tbp",
"cam_RR_A_0",
"cam_RR_A_rbp",
"cam_RR_A_tbp",
"cam_RR_B_0",
"cam_RR_B_rbp",
"cam_RR_B_tbp",
"cam_RR_C_0",
"cam_RR_C_rbp",
"cam_RR_C_tbp",
"cam_RR_D_0",
"cam_RR_D_rbp",
"cam_RR_D_tbp",
"cam_RR_AB_0",
"cam_RR_AB_rbp",
"cam_RR_AB_tbp",
"cam_RR_AC_0",
"cam_RR_AC_rbp",
"cam_RR_AC_tbp",
"cam_RR_AD_0",
"cam_RR_AD_rbp",
"cam_RR_AD_tbp",
"cam_RR_BC_0",
"cam_RR_BC_rbp",
"cam_RR_BC_tbp",
"cam_RR_BD_0",
"cam_RR_BD_rbp",
"cam_RR_BD_tbp",
"cam_RR_CD_0",
"cam_RR_CD_rbp",
"cam_RR_CD_tbp",
"cam_RR_ABC_0",
"cam_RR_ABC_rbp",
"cam_RR_ABC_tbp",
"cam_RR_ABD_0",
"cam_RR_ABD_rbp",
"cam_RR_ABD_tbp",
"cam_RR_ACD_0",
"cam_RR_ACD_rbp",
"cam_RR_ACD_tbp",
"cam_RR_BCD_0",
"cam_RR_BCD_rbp",
"cam_RR_BCD_tbp",
"cam_RR_ABCD_0",
"cam_RR_ABCD_rbp",
"cam_RR_ABCD_tbp",
"cam_RT_0_0",
"cam_RT_0_rbp",
"cam_RT_0_tbp",
"cam_RT_A_0",
"cam_RT_A_rbp",
"cam_RT_A_tbp",
"cam_RT_B_0",
"cam_RT_B_rbp",
"cam_RT_B_tbp",
"cam_RT_C_0",
"cam_RT_C_rbp",
"cam_RT_C_tbp",
"cam_RT_D_0",
"cam_RT_D_rbp",
"cam_RT_D_tbp",
"cam_RT_AB_0",
"cam_RT_AB_rbp",
"cam_RT_AB_tbp",
"cam_RT_AC_0",
"cam_RT_AC_rbp",
"cam_RT_AC_tbp",
"cam_RT_AD_0",
"cam_RT_AD_rbp",
"cam_RT_AD_tbp",
"cam_RT_BC_0",
"cam_RT_BC_rbp",
"cam_RT_BC_tbp",
"cam_RT_BD_0",
"cam_RT_BD_rbp",
"cam_RT_BD_tbp",
"cam_RT_CD_0",
"cam_RT_CD_rbp",
"cam_RT_CD_tbp",
"cam_RT_ABC_0",
"cam_RT_ABC_rbp",
"cam_RT_ABC_tbp",
"cam_RT_ABD_0",
"cam_RT_ABD_rbp",
"cam_RT_ABD_tbp",
"cam_RT_ACD_0",
"cam_RT_ACD_rbp",
"cam_RT_ACD_tbp",
"cam_RT_BCD_0",
"cam_RT_BCD_rbp",
"cam_RT_BCD_tbp",
"cam_RT_ABCD_0",
"cam_RT_ABCD_rbp",
"cam_RT_ABCD_tbp",
"cam_TR_0_0",
"cam_TR_0_rbp",
"cam_TR_0_tbp",
"cam_TR_A_0",
"cam_TR_A_rbp",
"cam_TR_A_tbp",
"cam_TR_B_0",
"cam_TR_B_rbp",
"cam_TR_B_tbp",
"cam_TR_C_0",
"cam_TR_C_rbp",
"cam_TR_C_tbp",
"cam_TR_D_0",
"cam_TR_D_rbp",
"cam_TR_D_tbp",
"cam_TR_AB_0",
"cam_TR_AB_rbp",
"cam_TR_AB_tbp",
"cam_TR_AC_0",
"cam_TR_AC_rbp",
"cam_TR_AC_tbp",
"cam_TR_AD_0",
"cam_TR_AD_rbp",
"cam_TR_AD_tbp",
"cam_TR_BC_0",
"cam_TR_BC_rbp",
"cam_TR_BC_tbp",
"cam_TR_BD_0",
"cam_TR_BD_rbp",
"cam_TR_BD_tbp",
"cam_TR_CD_0",
"cam_TR_CD_rbp",
"cam_TR_CD_tbp",
"cam_TR_ABC_0",
"cam_TR_ABC_rbp",
"cam_TR_ABC_tbp",
"cam_TR_ABD_0",
"cam_TR_ABD_rbp",
"cam_TR_ABD_tbp",
"cam_TR_ACD_0",
"cam_TR_ACD_rbp",
"cam_TR_ACD_tbp",
"cam_TR_BCD_0",
"cam_TR_BCD_rbp",
"cam_TR_BCD_tbp",
"cam_TR_ABCD_0",
"cam_TR_ABCD_rbp",
"cam_TR_ABCD_tbp",
"cam_TT_0_0",
"cam_TT_0_rbp",
"cam_TT_0_tbp",
"cam_TT_A_0",
"cam_TT_A_rbp",
"cam_TT_A_tbp",
"cam_TT_B_0",
"cam_TT_B_rbp",
"cam_TT_B_tbp",
"cam_TT_C_0",
"cam_TT_C_rbp",
"cam_TT_C_tbp",
"cam_TT_D_0",
"cam_TT_D_rbp",
"cam_TT_D_tbp",
"cam_TT_AB_0",
"cam_TT_AB_rbp",
"cam_TT_AB_tbp",
"cam_TT_AC_0",
"cam_TT_AC_rbp",
"cam_TT_AC_tbp",
"cam_TT_AD_0",
"cam_TT_AD_rbp",
"cam_TT_AD_tbp",
"cam_TT_BC_0",
"cam_TT_BC_rbp",
"cam_TT_BC_tbp",
"cam_TT_BD_0",
"cam_TT_BD_rbp",
"cam_TT_BD_tbp",
"cam_TT_CD_0",
"cam_TT_CD_rbp",
"cam_TT_CD_tbp",
"cam_TT_ABC_0",
"cam_TT_ABC_rbp",
"cam_TT_ABC_tbp",
"cam_TT_ABD_0",
"cam_TT_ABD_rbp",
"cam_TT_ABD_tbp",
"cam_TT_ACD_0",
"cam_TT_ACD_rbp",
"cam_TT_ACD_tbp",
"cam_TT_BCD_0",
"cam_TT_BCD_rbp",
"cam_TT_BCD_tbp",
"cam_TT_ABCD_0",
"cam_TT_ABCD_rbp",
"cam_TT_ABCD_tbp",
"ca",
"rbp",
"tbp",};
    }

    /**
     * @brief Get names of the fixed model parameters
     * @return the names
     */
    virtual std::vector<std::string> getFixedParameterNames() const override {
        return std::vector<std::string>{
            };
    }

    /**
     * @brief Get names of the observables
     * @return the names
     */
    virtual std::vector<std::string> getObservableNames() const override {
        return std::vector<std::string>{"x0",
"x1",
"x2",
"x3",
"x4",
"x5",
"x6",
"x7",
"x8",
"x9",
"x10",
"x11",
"x12",
"x13",
"x14",
"x15",
"x16",
"x17",
"x18",
"x19",
"x20",
"x21",
"x22",
"x23",
"x24",
"x25",
"x26",
"x27",
"x28",
"x29",
"x30",
"x31",
"x32",
"x33",
"x34",
"x35",
"x36",
"x37",
"x38",
"x39",
"x40",
"x41",
"x42",
"x43",
"x44",
"x45",
"x46",
"x47",
"x48",
"x49",
"x50",
"x51",
"x52",
"x53",
"x54",
"x55",
"x56",
"x57",
"x58",
"x59",
"x60",
"x61",
"x62",
"x63",
"x64",
"x65",
"x66",
"x67",
"x68",
"x69",
"x70",
"x71",
"x72",
"x73",
"x74",
"x75",
"x76",
"x77",
"x78",
"x79",
"x80",
"x81",
"x82",
"x83",
"x84",
"x85",
"x86",
"x87",
"x88",
"x89",
"x90",
"x91",
"x92",
"x93",
"x94",
"x95",
"x96",
"x97",
"x98",
"x99",
"x100",
"x101",
"x102",
"x103",
"x104",
"x105",
"x106",
"x107",
"x108",
"x109",
"x110",
"x111",
"x112",
"x113",
"x114",
"x115",
"x116",
"x117",
"x118",
"x119",
"x120",
"x121",
"x122",
"x123",
"x124",
"x125",
"x126",
"x127",
"x128",
"x129",
"x130",
"x131",
"x132",
"x133",
"x134",
"x135",
"x136",
"x137",
"x138",
"x139",
"x140",
"x141",
"x142",
"x143",
"x144",
"x145",
"x146",
"x147",
"x148",
"x149",
"x150",
"x151",
"x152",
"x153",
"x154",
"x155",
"x156",
"x157",
"x158",
"x159",
"x160",
"x161",
"x162",
"x163",
"x164",
"x165",
"x166",
"x167",
"x168",
"x169",
"x170",
"x171",
"x172",
"x173",
"x174",
"x175",
"x176",
"x177",
"x178",
"x179",
"x180",
"x181",
"x182",
"x183",
"x184",
"x185",
"x186",
"x187",
"x188",
"x189",
"x190",
"x191",
"x192",
"x193",
"x194",};
    }

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getParameterIds() const override {
        return std::vector<std::string>{"lC",
"lN",
"Kd_rbp_TT",
"kon_tbp",
"Kd_rbp_RT",
"kon_AT",
"kon_AR",
"kon_CR",
"cN",
"k_R2T_C",
"cC",
"kon_CT",
"k_R2T_N",
"Kd_tbp_RT",
"KDT",
"conc_rbp",
"KBT",
"Kd_rbp_RR",
"Kd_rbp_TR",
"Kd_tbp_TT",
"Kd_tbp_TR",
"conc_cam",
"Kd_tbp_RR",
"conc_tbp",
"kon_rbp",
"KCT",
"KAT",};
    }

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    virtual std::vector<std::string> getStateIds() const override {
        return std::vector<std::string>{"cam_RR_0_0",
"cam_RR_0_rbp",
"cam_RR_0_tbp",
"cam_RR_A_0",
"cam_RR_A_rbp",
"cam_RR_A_tbp",
"cam_RR_B_0",
"cam_RR_B_rbp",
"cam_RR_B_tbp",
"cam_RR_C_0",
"cam_RR_C_rbp",
"cam_RR_C_tbp",
"cam_RR_D_0",
"cam_RR_D_rbp",
"cam_RR_D_tbp",
"cam_RR_AB_0",
"cam_RR_AB_rbp",
"cam_RR_AB_tbp",
"cam_RR_AC_0",
"cam_RR_AC_rbp",
"cam_RR_AC_tbp",
"cam_RR_AD_0",
"cam_RR_AD_rbp",
"cam_RR_AD_tbp",
"cam_RR_BC_0",
"cam_RR_BC_rbp",
"cam_RR_BC_tbp",
"cam_RR_BD_0",
"cam_RR_BD_rbp",
"cam_RR_BD_tbp",
"cam_RR_CD_0",
"cam_RR_CD_rbp",
"cam_RR_CD_tbp",
"cam_RR_ABC_0",
"cam_RR_ABC_rbp",
"cam_RR_ABC_tbp",
"cam_RR_ABD_0",
"cam_RR_ABD_rbp",
"cam_RR_ABD_tbp",
"cam_RR_ACD_0",
"cam_RR_ACD_rbp",
"cam_RR_ACD_tbp",
"cam_RR_BCD_0",
"cam_RR_BCD_rbp",
"cam_RR_BCD_tbp",
"cam_RR_ABCD_0",
"cam_RR_ABCD_rbp",
"cam_RR_ABCD_tbp",
"cam_RT_0_0",
"cam_RT_0_rbp",
"cam_RT_0_tbp",
"cam_RT_A_0",
"cam_RT_A_rbp",
"cam_RT_A_tbp",
"cam_RT_B_0",
"cam_RT_B_rbp",
"cam_RT_B_tbp",
"cam_RT_C_0",
"cam_RT_C_rbp",
"cam_RT_C_tbp",
"cam_RT_D_0",
"cam_RT_D_rbp",
"cam_RT_D_tbp",
"cam_RT_AB_0",
"cam_RT_AB_rbp",
"cam_RT_AB_tbp",
"cam_RT_AC_0",
"cam_RT_AC_rbp",
"cam_RT_AC_tbp",
"cam_RT_AD_0",
"cam_RT_AD_rbp",
"cam_RT_AD_tbp",
"cam_RT_BC_0",
"cam_RT_BC_rbp",
"cam_RT_BC_tbp",
"cam_RT_BD_0",
"cam_RT_BD_rbp",
"cam_RT_BD_tbp",
"cam_RT_CD_0",
"cam_RT_CD_rbp",
"cam_RT_CD_tbp",
"cam_RT_ABC_0",
"cam_RT_ABC_rbp",
"cam_RT_ABC_tbp",
"cam_RT_ABD_0",
"cam_RT_ABD_rbp",
"cam_RT_ABD_tbp",
"cam_RT_ACD_0",
"cam_RT_ACD_rbp",
"cam_RT_ACD_tbp",
"cam_RT_BCD_0",
"cam_RT_BCD_rbp",
"cam_RT_BCD_tbp",
"cam_RT_ABCD_0",
"cam_RT_ABCD_rbp",
"cam_RT_ABCD_tbp",
"cam_TR_0_0",
"cam_TR_0_rbp",
"cam_TR_0_tbp",
"cam_TR_A_0",
"cam_TR_A_rbp",
"cam_TR_A_tbp",
"cam_TR_B_0",
"cam_TR_B_rbp",
"cam_TR_B_tbp",
"cam_TR_C_0",
"cam_TR_C_rbp",
"cam_TR_C_tbp",
"cam_TR_D_0",
"cam_TR_D_rbp",
"cam_TR_D_tbp",
"cam_TR_AB_0",
"cam_TR_AB_rbp",
"cam_TR_AB_tbp",
"cam_TR_AC_0",
"cam_TR_AC_rbp",
"cam_TR_AC_tbp",
"cam_TR_AD_0",
"cam_TR_AD_rbp",
"cam_TR_AD_tbp",
"cam_TR_BC_0",
"cam_TR_BC_rbp",
"cam_TR_BC_tbp",
"cam_TR_BD_0",
"cam_TR_BD_rbp",
"cam_TR_BD_tbp",
"cam_TR_CD_0",
"cam_TR_CD_rbp",
"cam_TR_CD_tbp",
"cam_TR_ABC_0",
"cam_TR_ABC_rbp",
"cam_TR_ABC_tbp",
"cam_TR_ABD_0",
"cam_TR_ABD_rbp",
"cam_TR_ABD_tbp",
"cam_TR_ACD_0",
"cam_TR_ACD_rbp",
"cam_TR_ACD_tbp",
"cam_TR_BCD_0",
"cam_TR_BCD_rbp",
"cam_TR_BCD_tbp",
"cam_TR_ABCD_0",
"cam_TR_ABCD_rbp",
"cam_TR_ABCD_tbp",
"cam_TT_0_0",
"cam_TT_0_rbp",
"cam_TT_0_tbp",
"cam_TT_A_0",
"cam_TT_A_rbp",
"cam_TT_A_tbp",
"cam_TT_B_0",
"cam_TT_B_rbp",
"cam_TT_B_tbp",
"cam_TT_C_0",
"cam_TT_C_rbp",
"cam_TT_C_tbp",
"cam_TT_D_0",
"cam_TT_D_rbp",
"cam_TT_D_tbp",
"cam_TT_AB_0",
"cam_TT_AB_rbp",
"cam_TT_AB_tbp",
"cam_TT_AC_0",
"cam_TT_AC_rbp",
"cam_TT_AC_tbp",
"cam_TT_AD_0",
"cam_TT_AD_rbp",
"cam_TT_AD_tbp",
"cam_TT_BC_0",
"cam_TT_BC_rbp",
"cam_TT_BC_tbp",
"cam_TT_BD_0",
"cam_TT_BD_rbp",
"cam_TT_BD_tbp",
"cam_TT_CD_0",
"cam_TT_CD_rbp",
"cam_TT_CD_tbp",
"cam_TT_ABC_0",
"cam_TT_ABC_rbp",
"cam_TT_ABC_tbp",
"cam_TT_ABD_0",
"cam_TT_ABD_rbp",
"cam_TT_ABD_tbp",
"cam_TT_ACD_0",
"cam_TT_ACD_rbp",
"cam_TT_ACD_tbp",
"cam_TT_BCD_0",
"cam_TT_BCD_rbp",
"cam_TT_BCD_tbp",
"cam_TT_ABCD_0",
"cam_TT_ABCD_rbp",
"cam_TT_ABCD_tbp",
"ca",
"rbp",
"tbp",};
    }

    /**
     * @brief Get ids of the fixed model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getFixedParameterIds() const override {
        return std::vector<std::string>{
            };
    }

    /**
     * @brief Get ids of the observables
     * @return the ids
     */
    virtual std::vector<std::string> getObservableIds() const override {
        return std::vector<std::string>{"y0",
"y1",
"y2",
"y3",
"y4",
"y5",
"y6",
"y7",
"y8",
"y9",
"y10",
"y11",
"y12",
"y13",
"y14",
"y15",
"y16",
"y17",
"y18",
"y19",
"y20",
"y21",
"y22",
"y23",
"y24",
"y25",
"y26",
"y27",
"y28",
"y29",
"y30",
"y31",
"y32",
"y33",
"y34",
"y35",
"y36",
"y37",
"y38",
"y39",
"y40",
"y41",
"y42",
"y43",
"y44",
"y45",
"y46",
"y47",
"y48",
"y49",
"y50",
"y51",
"y52",
"y53",
"y54",
"y55",
"y56",
"y57",
"y58",
"y59",
"y60",
"y61",
"y62",
"y63",
"y64",
"y65",
"y66",
"y67",
"y68",
"y69",
"y70",
"y71",
"y72",
"y73",
"y74",
"y75",
"y76",
"y77",
"y78",
"y79",
"y80",
"y81",
"y82",
"y83",
"y84",
"y85",
"y86",
"y87",
"y88",
"y89",
"y90",
"y91",
"y92",
"y93",
"y94",
"y95",
"y96",
"y97",
"y98",
"y99",
"y100",
"y101",
"y102",
"y103",
"y104",
"y105",
"y106",
"y107",
"y108",
"y109",
"y110",
"y111",
"y112",
"y113",
"y114",
"y115",
"y116",
"y117",
"y118",
"y119",
"y120",
"y121",
"y122",
"y123",
"y124",
"y125",
"y126",
"y127",
"y128",
"y129",
"y130",
"y131",
"y132",
"y133",
"y134",
"y135",
"y136",
"y137",
"y138",
"y139",
"y140",
"y141",
"y142",
"y143",
"y144",
"y145",
"y146",
"y147",
"y148",
"y149",
"y150",
"y151",
"y152",
"y153",
"y154",
"y155",
"y156",
"y157",
"y158",
"y159",
"y160",
"y161",
"y162",
"y163",
"y164",
"y165",
"y166",
"y167",
"y168",
"y169",
"y170",
"y171",
"y172",
"y173",
"y174",
"y175",
"y176",
"y177",
"y178",
"y179",
"y180",
"y181",
"y182",
"y183",
"y184",
"y185",
"y186",
"y187",
"y188",
"y189",
"y190",
"y191",
"y192",
"y193",
"y194",};
    }

    /**
     * @brief function indicating whether reinitialization of states depending on
     fixed parameters is permissible
     * @return flag inidication whether reinitialization of states depending on
     fixed parameters is permissible
     */
    virtual bool isFixedParameterStateReinitializationAllowed() const override {
        return true;
    }

    /**
     * @brief returns the amici version that was used to generate the model
     * @return ver amici version string
     */
    virtual const std::string getAmiciVersion() const override {
        return "0.10.7";
    }

    /**
     & @brief returns the amici version that was used to generate the model
     * @return commit amici git commit hash
     */
    virtual const std::string getAmiciCommit() const override {
        return "unknown";
    }

    virtual bool wasPythonGenerated() const override {
        return true;
    }
};

#endif /* _amici_TPL_MODELNAME_h */
