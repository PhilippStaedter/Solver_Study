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
extern void J_Pathak2013a(realtype *J, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w,
                            const realtype *dwdx);
extern void JB_Pathak2013a(realtype *JB, const realtype t, const realtype *x,
                             const realtype *p, const realtype *k,
                             const realtype *h, const realtype *xB,
                             const realtype *w, const realtype *dwdx);
extern void JDiag_Pathak2013a(realtype *JDiag, const realtype t,
                                const realtype *x, const realtype *p,
                                const realtype *k, const realtype *h,
                                const realtype *w, const realtype *dwdx);
extern void JSparse_Pathak2013a(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JSparse_colptrs_Pathak2013a(sunindextype *colptrs);
extern void JSparse_rowvals_Pathak2013a(sunindextype *rowvals);
extern void JSparseB_Pathak2013a(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void JSparseB_colptrs_Pathak2013a(sunindextype *colptrs);
extern void JSparseB_rowvals_Pathak2013a(sunindextype *rowvals);
extern void Jy_Pathak2013a(realtype *nllh, const int iy, const realtype *p,
                             const realtype *k, const realtype *y,
                             const realtype *sigmay, const realtype *my);
extern void dJydsigmay_Pathak2013a(realtype *dJydsigmay, const int iy,
                                     const realtype *p, const realtype *k,
                                     const realtype *y, const realtype *sigmay,
                                     const realtype *my);
extern void dJydy_Pathak2013a(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_Pathak2013a(sunindextype *colptrs, int index);
extern void dJydy_rowvals_Pathak2013a(sunindextype *rowvals, int index);
extern void dwdp_Pathak2013a(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip);
extern void dwdx_Pathak2013a(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl);
extern void dwdx_colptrs_Pathak2013a(sunindextype *colptrs);
extern void dwdx_rowvals_Pathak2013a(sunindextype *rowvals);
extern void dxdotdw_Pathak2013a(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdw_colptrs_Pathak2013a(sunindextype *colptrs);
extern void dxdotdw_rowvals_Pathak2013a(sunindextype *rowvals);
extern void dxdotdp_Pathak2013a(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w);
extern void dydx_Pathak2013a(realtype *dydx, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w, const realtype *dwdx);
extern void dydp_Pathak2013a(realtype *dydp, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const int ip, const realtype *w,
                               const realtype *dwp);
extern void dsigmaydp_Pathak2013a(realtype *dsigmaydp, const realtype t,
                                    const realtype *p, const realtype *k,
                                    const int ip);
extern void sigmay_Pathak2013a(realtype *sigmay, const realtype t,
                                 const realtype *p, const realtype *k);
extern void w_Pathak2013a(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl);
extern void x0_Pathak2013a(realtype *x0, const realtype t, const realtype *p,
                             const realtype *k);
extern void x0_fixedParameters_Pathak2013a(realtype *x0, const realtype t,
                                             const realtype *p,
                                             const realtype *k);
extern void sx0_Pathak2013a(realtype *sx0, const realtype t,
                              const realtype *x0, const realtype *p,
                              const realtype *k, const int ip);
extern void sx0_fixedParameters_Pathak2013a(realtype *sx0, const realtype t,
                                              const realtype *x0,
                                              const realtype *p,
                                              const realtype *k, const int ip);
extern void xdot_Pathak2013a(realtype *xdot, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w);
extern void y_Pathak2013a(realtype *y, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w);

extern void x_solver_Pathak2013a(realtype *x_solver, const realtype *x_rdata);
extern void total_cl_Pathak2013a(realtype *total_cl, const realtype *x_rdata);

/**
 * @brief AMICI-generated model subclass.
 */
class Model_Pathak2013a : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_Pathak2013a()
        : amici::Model_ODE(
              52,                                // nx_rdata
              52,                            // nxtrue_rdata
              52,                               // nx_solver
              52,                           // nxtrue_solver
              52,                                      // ny
              52,                                  // nytrue
              0,                                      // nz
              0,                                  // nztrue
              0,                                  // nevent
              1,                              // nobjective
              88,                                      // nw
              176,                                   // ndwdx
              15488,                                   // ndwdp
              176,                                // ndxdotdw
              std::vector<int>{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},                                  // ndjydy
              228,                                     // nnz
              52,                                     // ubw
              52,                                     // lbw
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},       // dynamic parameters
              std::vector<realtype>{}, // fixedParameters
              std::vector<int>{},                          // plist
              std::vector<realtype>(52, 0.0),   // idlist
              std::vector<int>{}                           // z2event
          ) {}

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    virtual amici::Model *clone() const override {
        return new Model_Pathak2013a(*this);
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
        J_Pathak2013a(J, t, x, p, k, h, w, dwdx);
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
        JB_Pathak2013a(JB, t, x, p, k, h, xB, w, dwdx);
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
        JDiag_Pathak2013a(JDiag, t, x, p, k, h, w, dwdx);
    }

    virtual void fJSparse(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        JSparse_Pathak2013a( JSparse,  t,  x,  p,  k,  h,  w,  dwdx);
    }


    virtual void fJSparse_colptrs(sunindextype *colptrs) override {
        JSparse_colptrs_Pathak2013a(colptrs);
    }


    virtual void fJSparse_rowvals(sunindextype *rowvals) override {
        JSparse_rowvals_Pathak2013a(rowvals);
    }


    virtual void fJSparseB(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        JSparseB_Pathak2013a( JSparseB,  t,  x,  p,  k,  h,  xB,  w,  dwdx);
    }


    virtual void fJSparseB_colptrs(sunindextype *colptrs) override {
        JSparseB_colptrs_Pathak2013a(colptrs);
    }


    virtual void fJSparseB_rowvals(sunindextype *rowvals) override {
        JSparseB_rowvals_Pathak2013a(rowvals);
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
        Jy_Pathak2013a(nllh, iy, p, k, y, sigmay, my);
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
        dJydsigmay_Pathak2013a(dJydsigma, iy, p, k, y, sigmay, my);
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
        dsigmaydp_Pathak2013a(dsigmaydp, t, p, k, ip);
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
        dJydy_Pathak2013a( dJydy,  iy,  p,  k,  y,  sigmay,  my);
    }

    virtual void fdJydy_colptrs(sunindextype *colptrs, int index) override {
        dJydy_colptrs_Pathak2013a(colptrs, index);
    }

    virtual void fdJydy_rowvals(sunindextype *rowvals, int index) override {
        dJydy_rowvals_Pathak2013a(rowvals, index);
    }


    virtual void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip) override {
        dwdp_Pathak2013a( dwdp,  t,  x,  p,  k,  h,  w,  tcl,  dtcldp,  ip);
    }


    virtual void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {
        dwdx_Pathak2013a( dwdx,  t,  x,  p,  k,  h,  w,  tcl);
    }


    virtual void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdw_Pathak2013a( dxdotdw,  t,  x,  p,  k,  h,  w);
    }


    virtual void fdxdotdw_colptrs(sunindextype *colptrs) override {
        dxdotdw_colptrs_Pathak2013a(colptrs);
    }


    virtual void fdxdotdw_rowvals(sunindextype *rowvals) override {
        dxdotdw_rowvals_Pathak2013a(rowvals);
    }


    virtual void fdxdotdp(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w) override {
        dxdotdp_Pathak2013a( dxdotdp,  t,  x,  p,  k,  h,  ip,  w);
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
        dydx_Pathak2013a(dydx, t, x, p, k, h, w, dwdx);
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
        dydp_Pathak2013a(dydp, t, x, p, k, h, ip, w, dwdp);
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
        sigmay_Pathak2013a(sigmay, t, p, k);
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
        sx0_Pathak2013a(sx0, t, x0, p, k, ip);
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
        sx0_fixedParameters_Pathak2013a(sx0, t, x0, p, k, ip);
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
        w_Pathak2013a( w,  t,  x,  p,  k,  h,  tcl);
    }


    /** model specific implementation of fx0
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0(realtype *x0, const realtype t, const realtype *p,
                     const realtype *k) override {
        x0_Pathak2013a(x0, t, p, k);
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
        x0_fixedParameters_Pathak2013a(x0, t, p, k);
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
        xdot_Pathak2013a(xdot, t, x, p, k, h, w);
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
        y_Pathak2013a(y, t, x, p, k, h, w);
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
        x_solver_Pathak2013a( x_solver,  x_rdata);
    }


    virtual void ftotal_cl(realtype *total_cl, const realtype *x_rdata) override {
        total_cl_Pathak2013a( total_cl,  x_rdata);
    }


    /**
     * @brief Get names of the model parameters
     * @return the names
     */
    virtual std::vector<std::string> getParameterNames() const override {
        return std::vector<std::string>{"Dissociation constant of reaction re1",
"Association constant of reaction re1",
"Dissociation constant of reaction re2",
"Association constant of reaction re2",
"Dissociation constant of reaction re3",
"Association constant of reaction re3",
"Dissociation constant of reaction re4",
"Association constant of reaction re4",
"Dissociation constant of reaction re5",
"Association constant of reaction re5",
"Dissociation constant of reaction re6",
"Association constant of reaction re6",
"Dissociation constant of reaction re7",
"Association constant of reaction re7",
"Dissociation constant of reaction re8",
"Association constant of reaction re8",
"Dissociation constant of reaction re9",
"Association constant of reaction re9",
"Dissociation constant of reaction re10",
"Association constant of reaction re10",
"Dissociation constant of reaction re11",
"Association constant of reaction re11",
"Dissociation constant of reaction re12",
"Association constant of reaction re12",
"Dissociation constant of reaction re13",
"Association constant of reaction re13",
"Dissociation constant of reaction re14",
"Association constant of reaction re14",
"Dissociation constant of reaction re15",
"Association constant of reaction re15",
"Dissociation constant of reaction re16",
"Association constant of reaction re16",
"Dissociation constant of reaction re17",
"Association constant of reaction re17",
"Dissociation constant of reaction re18",
"Association constant of reaction re18",
"Dissociation constant of reaction re19",
"Association constant of reaction re19",
"Dissociation constant of reaction re20",
"Association constant of reaction re20",
"Dissociation constant of reaction re21",
"Association constant of reaction re21",
"Dissociation constant of reaction re22",
"Association constant of reaction re22",
"Dissociation constant of reaction re23",
"Association constant of reaction re23",
"Dissociation constant of reaction re24",
"Association constant of reaction re24",
"Dissociation constant of reaction re25",
"Association constant of reaction re25",
"Dissociation constant of reaction re26",
"Association constant of reaction re26",
"Dissociation constant of reaction re27",
"Association constant of reaction re27",
"Dissociation constant of reaction re28",
"Association constant of reaction re28",
"Dissociation constant of reaction re29",
"Association constant of reaction re29",
"Dissociation constant of reaction re30",
"Association constant of reaction re30",
"Dissociation constant of reaction re31",
"Association constant of reaction re31",
"Dissociation constant of reaction re32",
"Association constant of reaction re32",
"Dissociation constant of reaction re33",
"Association constant of reaction re33",
"Dissociation constant of reaction re34",
"Association constant of reaction re34",
"Dissociation constant of reaction re35",
"Association constant of reaction re35",
"Dissociation constant of reaction re36",
"Association constant of reaction re36",
"Dissociation constant of reaction re37",
"Association constant of reaction re37",
"Dissociation constant of reaction re38",
"Association constant of reaction re38",
"Dissociation constant of reaction re39",
"Association constant of reaction re39",
"Dissociation constant of reaction re40",
"Association constant of reaction re40",
"Dissociation constant of reaction re41",
"Association constant of reaction re41",
"Dissociation constant of reaction re42",
"Association constant of reaction re42",
"Dissociation constant of reaction re43",
"Association constant of reaction re43",
"Dissociation constant of reaction re44",
"Association constant of reaction re44",
"Dissociation constant of reaction re45",
"Association constant of reaction re45",
"Dissociation constant of reaction re46",
"Association constant of reaction re46",
"Dissociation constant of reaction re47",
"Association constant of reaction re47",
"Dissociation constant of reaction re48",
"Association constant of reaction re48",
"Dissociation constant of reaction re49",
"Association constant of reaction re49",
"Dissociation constant of reaction re50",
"Association constant of reaction re50",
"Dissociation constant of reaction re51",
"Association constant of reaction re51",
"Dissociation constant of reaction re52",
"Association constant of reaction re52",
"Dissociation constant of reaction re53",
"Association constant of reaction re53",
"Dissociation constant of reaction re54",
"Association constant of reaction re54",
"Dissociation constant of reaction re55",
"Association constant of reaction re55",
"Dissociation constant of reaction re56",
"Association constant of reaction re56",
"Dissociation constant of reaction re57",
"Association constant of reaction re57",
"Dissociation constant of reaction re58",
"Association constant of reaction re58",
"Dissociation constant of reaction re59",
"Association constant of reaction re59",
"Dissociation constant of reaction re60",
"Association constant of reaction re60",
"Dissociation constant of reaction re61",
"Association constant of reaction re61",
"Dissociation constant of reaction re62",
"Association constant of reaction re62",
"Dissociation constant of reaction re63",
"Association constant of reaction re63",
"Dissociation constant of reaction re64",
"Association constant of reaction re64",
"Dissociation constant of reaction re65",
"Association constant of reaction re65",
"Dissociation constant of reaction re66",
"Association constant of reaction re66",
"Dissociation constant of reaction re67",
"Association constant of reaction re67",
"Dissociation constant of reaction re68",
"Association constant of reaction re68",
"Dissociation constant of reaction re69",
"Association constant of reaction re69",
"Dissociation constant of reaction re70",
"Association constant of reaction re70",
"Dissociation constant of reaction re71",
"Association constant of reaction re71",
"Dissociation constant of reaction re72",
"Association constant of reaction re72",
"Dissociation constant of reaction re73",
"Association constant of reaction re73",
"Dissociation constant of reaction re74",
"Association constant of reaction re74",
"Dissociation constant of reaction re75",
"Association constant of reaction re75",
"Dissociation constant of reaction re76",
"Association constant of reaction re76",
"Dissociation constant of reaction re77",
"Association constant of reaction re77",
"Dissociation constant of reaction re78",
"Association constant of reaction re78",
"Dissociation constant of reaction re79",
"Association constant of reaction re79",
"Dissociation constant of reaction re80",
"Association constant of reaction re80",
"Dissociation constant of reaction re81",
"Association constant of reaction re81",
"Dissociation constant of reaction re82",
"Association constant of reaction re82",
"Dissociation constant of reaction re83",
"Association constant of reaction re83",
"Dissociation constant of reaction re84",
"Association constant of reaction re84",
"Dissociation constant of reaction re85",
"Association constant of reaction re85",
"Dissociation constant of reaction re86",
"Association constant of reaction re86",
"Dissociation constant of reaction re87",
"Association constant of reaction re87",
"Dissociation constant of reaction re88",
"Association constant of reaction re88",};
    }

    /**
     * @brief Get names of the model states
     * @return the names
     */
    virtual std::vector<std::string> getStateNames() const override {
        return std::vector<std::string>{"Fungal pathogen",
"Bacterial pathogen",
"LysM",
"PRRs",
"FLS2",
"LRR",
"MAPKKK",
"MAPKKK",
"MAPKKK1",
"MAPKKK18",
"MAPKKK19",
"MAPKKK20",
"EDR1",
"MAPKK",
"MAPKK",
"MAPKK2",
"MAPKK4",
"MAPKK5",
"MAPKK9",
"MAPK",
"MAPK",
"MAPK2",
"MAPK3",
"MAPK4",
"MAPK6",
"WRKY1",
"WRKY1",
"MYB2",
"MYB2",
"WRKY33",
"WRKY33",
"WRKY6",
"WRKY6",
"MYB4",
"MYB4",
"WRKY25",
"WRKY25",
"WRKY12",
"WRKY12",
"WRKY22",
"WRKY22",
"WRKY28",
"WRKY28",
"WRKY29",
"WRKY29",
"MYB44",
"NAC",
"bZIP",
"AP2",
"Response",
"SIMK",
"SAMK",};
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
"x51",};
    }

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getParameterIds() const override {
        return std::vector<std::string>{"kdiss_re1",
"kass_re1",
"kdiss_re2",
"kass_re2",
"kdiss_re3",
"kass_re3",
"kdiss_re4",
"kass_re4",
"kdiss_re5",
"kass_re5",
"kdiss_re6",
"kass_re6",
"kdiss_re7",
"kass_re7",
"kdiss_re8",
"kass_re8",
"kdiss_re9",
"kass_re9",
"kdiss_re10",
"kass_re10",
"kdiss_re11",
"kass_re11",
"kdiss_re12",
"kass_re12",
"kdiss_re13",
"kass_re13",
"kdiss_re14",
"kass_re14",
"kdiss_re15",
"kass_re15",
"kdiss_re16",
"kass_re16",
"kdiss_re17",
"kass_re17",
"kdiss_re18",
"kass_re18",
"kdiss_re19",
"kass_re19",
"kdiss_re20",
"kass_re20",
"kdiss_re21",
"kass_re21",
"kdiss_re22",
"kass_re22",
"kdiss_re23",
"kass_re23",
"kdiss_re24",
"kass_re24",
"kdiss_re25",
"kass_re25",
"kdiss_re26",
"kass_re26",
"kdiss_re27",
"kass_re27",
"kdiss_re28",
"kass_re28",
"kdiss_re29",
"kass_re29",
"kdiss_re30",
"kass_re30",
"kdiss_re31",
"kass_re31",
"kdiss_re32",
"kass_re32",
"kdiss_re33",
"kass_re33",
"kdiss_re34",
"kass_re34",
"kdiss_re35",
"kass_re35",
"kdiss_re36",
"kass_re36",
"kdiss_re37",
"kass_re37",
"kdiss_re38",
"kass_re38",
"kdiss_re39",
"kass_re39",
"kdiss_re40",
"kass_re40",
"kdiss_re41",
"kass_re41",
"kdiss_re42",
"kass_re42",
"kdiss_re43",
"kass_re43",
"kdiss_re44",
"kass_re44",
"kdiss_re45",
"kass_re45",
"kdiss_re46",
"kass_re46",
"kdiss_re47",
"kass_re47",
"kdiss_re48",
"kass_re48",
"kdiss_re49",
"kass_re49",
"kdiss_re50",
"kass_re50",
"kdiss_re51",
"kass_re51",
"kdiss_re52",
"kass_re52",
"kdiss_re53",
"kass_re53",
"kdiss_re54",
"kass_re54",
"kdiss_re55",
"kass_re55",
"kdiss_re56",
"kass_re56",
"kdiss_re57",
"kass_re57",
"kdiss_re58",
"kass_re58",
"kdiss_re59",
"kass_re59",
"kdiss_re60",
"kass_re60",
"kdiss_re61",
"kass_re61",
"kdiss_re62",
"kass_re62",
"kdiss_re63",
"kass_re63",
"kdiss_re64",
"kass_re64",
"kdiss_re65",
"kass_re65",
"kdiss_re66",
"kass_re66",
"kdiss_re67",
"kass_re67",
"kdiss_re68",
"kass_re68",
"kdiss_re69",
"kass_re69",
"kdiss_re70",
"kass_re70",
"kdiss_re71",
"kass_re71",
"kdiss_re72",
"kass_re72",
"kdiss_re73",
"kass_re73",
"kdiss_re74",
"kass_re74",
"kdiss_re75",
"kass_re75",
"kdiss_re76",
"kass_re76",
"kdiss_re77",
"kass_re77",
"kdiss_re78",
"kass_re78",
"kdiss_re79",
"kass_re79",
"kdiss_re80",
"kass_re80",
"kdiss_re81",
"kass_re81",
"kdiss_re82",
"kass_re82",
"kdiss_re83",
"kass_re83",
"kdiss_re84",
"kass_re84",
"kdiss_re85",
"kass_re85",
"kdiss_re86",
"kass_re86",
"kdiss_re87",
"kass_re87",
"kdiss_re88",
"kass_re88",};
    }

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    virtual std::vector<std::string> getStateIds() const override {
        return std::vector<std::string>{"s1",
"s2",
"s3",
"s4",
"s5",
"s6",
"s7",
"s8",
"s9",
"s10",
"s11",
"s12",
"s13",
"s14",
"s15",
"s16",
"s17",
"s18",
"s19",
"s20",
"s21",
"s22",
"s23",
"s24",
"s25",
"s28",
"s29",
"s30",
"s31",
"s32",
"s33",
"s34",
"s35",
"s36",
"s37",
"s38",
"s39",
"s40",
"s41",
"s42",
"s43",
"s44",
"s45",
"s46",
"s47",
"s48",
"s49",
"s50",
"s51",
"s52",
"s26",
"s27",};
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
"y51",};
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
