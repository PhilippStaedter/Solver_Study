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
extern void J_Eungdamrong2007(realtype *J, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w,
                            const realtype *dwdx);
extern void JB_Eungdamrong2007(realtype *JB, const realtype t, const realtype *x,
                             const realtype *p, const realtype *k,
                             const realtype *h, const realtype *xB,
                             const realtype *w, const realtype *dwdx);
extern void JDiag_Eungdamrong2007(realtype *JDiag, const realtype t,
                                const realtype *x, const realtype *p,
                                const realtype *k, const realtype *h,
                                const realtype *w, const realtype *dwdx);
extern void JSparse_Eungdamrong2007(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JSparse_colptrs_Eungdamrong2007(sunindextype *colptrs);
extern void JSparse_rowvals_Eungdamrong2007(sunindextype *rowvals);
extern void JSparseB_Eungdamrong2007(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void JSparseB_colptrs_Eungdamrong2007(sunindextype *colptrs);
extern void JSparseB_rowvals_Eungdamrong2007(sunindextype *rowvals);
extern void Jy_Eungdamrong2007(realtype *nllh, const int iy, const realtype *p,
                             const realtype *k, const realtype *y,
                             const realtype *sigmay, const realtype *my);
extern void dJydsigmay_Eungdamrong2007(realtype *dJydsigmay, const int iy,
                                     const realtype *p, const realtype *k,
                                     const realtype *y, const realtype *sigmay,
                                     const realtype *my);
extern void dJydy_Eungdamrong2007(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_Eungdamrong2007(sunindextype *colptrs, int index);
extern void dJydy_rowvals_Eungdamrong2007(sunindextype *rowvals, int index);
extern void dwdp_Eungdamrong2007(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip);
extern void dwdx_Eungdamrong2007(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl);
extern void dwdx_colptrs_Eungdamrong2007(sunindextype *colptrs);
extern void dwdx_rowvals_Eungdamrong2007(sunindextype *rowvals);
extern void dxdotdw_Eungdamrong2007(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdw_colptrs_Eungdamrong2007(sunindextype *colptrs);
extern void dxdotdw_rowvals_Eungdamrong2007(sunindextype *rowvals);
extern void dxdotdp_Eungdamrong2007(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w);
extern void dydx_Eungdamrong2007(realtype *dydx, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w, const realtype *dwdx);
extern void dydp_Eungdamrong2007(realtype *dydp, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const int ip, const realtype *w,
                               const realtype *dwp);
extern void dsigmaydp_Eungdamrong2007(realtype *dsigmaydp, const realtype t,
                                    const realtype *p, const realtype *k,
                                    const int ip);
extern void sigmay_Eungdamrong2007(realtype *sigmay, const realtype t,
                                 const realtype *p, const realtype *k);
extern void w_Eungdamrong2007(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl);
extern void x0_Eungdamrong2007(realtype *x0, const realtype t, const realtype *p,
                             const realtype *k);
extern void x0_fixedParameters_Eungdamrong2007(realtype *x0, const realtype t,
                                             const realtype *p,
                                             const realtype *k);
extern void sx0_Eungdamrong2007(realtype *sx0, const realtype t,
                              const realtype *x0, const realtype *p,
                              const realtype *k, const int ip);
extern void sx0_fixedParameters_Eungdamrong2007(realtype *sx0, const realtype t,
                                              const realtype *x0,
                                              const realtype *p,
                                              const realtype *k, const int ip);
extern void xdot_Eungdamrong2007(realtype *xdot, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w);
extern void y_Eungdamrong2007(realtype *y, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w);

extern void x_solver_Eungdamrong2007(realtype *x_solver, const realtype *x_rdata);
extern void total_cl_Eungdamrong2007(realtype *total_cl, const realtype *x_rdata);

/**
 * @brief AMICI-generated model subclass.
 */
class Model_Eungdamrong2007 : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_Eungdamrong2007()
        : amici::Model_ODE(
              46,                                // nx_rdata
              46,                            // nxtrue_rdata
              46,                               // nx_solver
              46,                           // nxtrue_solver
              46,                                      // ny
              46,                                  // nytrue
              0,                                      // nz
              0,                                  // nztrue
              0,                                  // nevent
              1,                              // nobjective
              43,                                      // nw
              101,                                   // ndwdx
              5117,                                   // ndwdp
              94,                                // ndxdotdw
              std::vector<int>{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},                                  // ndjydy
              174,                                     // nnz
              46,                                     // ubw
              46,                                     // lbw
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{0.00166112956810631, 0.01, 0.05, 10.0, 1072.0, 0.045, 2857.0, 0.0015, 0.01, 0.05, 10.0, 2.1, 0.11, 0.2, 0.0, 1032.0, 0.0, 0.5, 0.1, 0.0, 0.5, 0.01, 120.0, 0.0, 0.5, 0.1, 600.0, 0.0, 0.01, 120.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0001, 0.0, 0.1, 120.0, 0.0, 0.1, 0.5, 0.0, 1e-05, 0.015, 0.0, 1200.0, 0.0, 0.01, 120.0, 0.0, 0.0168, 0.025, 0.01, 120.0, 0.0, 1e-05, 0.000158489319246111, 0.0, 0.0, 0.0001, 0.0, 0.0001, 1.18, 0.0, 1200.0, 0.0, 1.188, 0.0, 0.1, 90.0, 0.0, 0.0, 0.0, 0.002, 0.0, 10.0, 50.0, 1.0, 600.0, 0.0, 0.0, 0.5, 0.0, 0.01, 1.0, 0.0, 1.0, 3.0, 1.0, 15.0, 0.0, 0.0, 0.0001, 0.0, 5.0, 10.0, 0.0, 0.0, 410.0, 0.0, 1000.0, 0.0, 0.0, 0.25, 0.0, 600.0, 0.0, 0.27, 0.0664, 0.0, 4.69, 0.8, 0.0, 3.16e-05, 0.0, 0.0, 0.0001, 0.0},       // dynamic parameters
              std::vector<realtype>{}, // fixedParameters
              std::vector<int>{},                          // plist
              std::vector<realtype>(46, 0.0),   // idlist
              std::vector<int>{}                           // z2event
          ) {}

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    virtual amici::Model *clone() const override {
        return new Model_Eungdamrong2007(*this);
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
        J_Eungdamrong2007(J, t, x, p, k, h, w, dwdx);
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
        JB_Eungdamrong2007(JB, t, x, p, k, h, xB, w, dwdx);
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
        JDiag_Eungdamrong2007(JDiag, t, x, p, k, h, w, dwdx);
    }

    virtual void fJSparse(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        JSparse_Eungdamrong2007( JSparse,  t,  x,  p,  k,  h,  w,  dwdx);
    }


    virtual void fJSparse_colptrs(sunindextype *colptrs) override {
        JSparse_colptrs_Eungdamrong2007(colptrs);
    }


    virtual void fJSparse_rowvals(sunindextype *rowvals) override {
        JSparse_rowvals_Eungdamrong2007(rowvals);
    }


    virtual void fJSparseB(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        JSparseB_Eungdamrong2007( JSparseB,  t,  x,  p,  k,  h,  xB,  w,  dwdx);
    }


    virtual void fJSparseB_colptrs(sunindextype *colptrs) override {
        JSparseB_colptrs_Eungdamrong2007(colptrs);
    }


    virtual void fJSparseB_rowvals(sunindextype *rowvals) override {
        JSparseB_rowvals_Eungdamrong2007(rowvals);
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
        Jy_Eungdamrong2007(nllh, iy, p, k, y, sigmay, my);
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
        dJydsigmay_Eungdamrong2007(dJydsigma, iy, p, k, y, sigmay, my);
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
        dsigmaydp_Eungdamrong2007(dsigmaydp, t, p, k, ip);
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
        dJydy_Eungdamrong2007( dJydy,  iy,  p,  k,  y,  sigmay,  my);
    }

    virtual void fdJydy_colptrs(sunindextype *colptrs, int index) override {
        dJydy_colptrs_Eungdamrong2007(colptrs, index);
    }

    virtual void fdJydy_rowvals(sunindextype *rowvals, int index) override {
        dJydy_rowvals_Eungdamrong2007(rowvals, index);
    }


    virtual void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip) override {
        dwdp_Eungdamrong2007( dwdp,  t,  x,  p,  k,  h,  w,  tcl,  dtcldp,  ip);
    }


    virtual void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {
        dwdx_Eungdamrong2007( dwdx,  t,  x,  p,  k,  h,  w,  tcl);
    }


    virtual void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdw_Eungdamrong2007( dxdotdw,  t,  x,  p,  k,  h,  w);
    }


    virtual void fdxdotdw_colptrs(sunindextype *colptrs) override {
        dxdotdw_colptrs_Eungdamrong2007(colptrs);
    }


    virtual void fdxdotdw_rowvals(sunindextype *rowvals) override {
        dxdotdw_rowvals_Eungdamrong2007(rowvals);
    }


    virtual void fdxdotdp(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w) override {
        dxdotdp_Eungdamrong2007( dxdotdp,  t,  x,  p,  k,  h,  ip,  w);
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
        dydx_Eungdamrong2007(dydx, t, x, p, k, h, w, dwdx);
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
        dydp_Eungdamrong2007(dydp, t, x, p, k, h, ip, w, dwdp);
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
        sigmay_Eungdamrong2007(sigmay, t, p, k);
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
        sx0_Eungdamrong2007(sx0, t, x0, p, k, ip);
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
        sx0_fixedParameters_Eungdamrong2007(sx0, t, x0, p, k, ip);
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
        w_Eungdamrong2007( w,  t,  x,  p,  k,  h,  tcl);
    }


    /** model specific implementation of fx0
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0(realtype *x0, const realtype t, const realtype *p,
                     const realtype *k) override {
        x0_Eungdamrong2007(x0, t, p, k);
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
        x0_fixedParameters_Eungdamrong2007(x0, t, p, k);
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
        xdot_Eungdamrong2007(xdot, t, x, p, k, h, w);
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
        y_Eungdamrong2007(y, t, x, p, k, h, w);
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
        x_solver_Eungdamrong2007( x_solver,  x_rdata);
    }


    virtual void ftotal_cl(realtype *total_cl, const realtype *x_rdata) override {
        total_cl_Eungdamrong2007( total_cl,  x_rdata);
    }


    /**
     * @brief Get names of the model parameters
     * @return the names
     */
    virtual std::vector<std::string> getParameterNames() const override {
        return std::vector<std::string>{"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",};
    }

    /**
     * @brief Get names of the model states
     * @return the names
     */
    virtual std::vector<std::string> getStateNames() const override {
        return std::vector<std::string>{"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",};
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
"x45",};
    }

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getParameterIds() const override {
        return std::vector<std::string>{"KMOLE",
"kStimSynPIP2_PIP2_synthesis",
"tauPIP2syn_PIP2_synthesis",
"PIP2syndecay_PIP2_synthesis",
"PIP2_basal_PIP2_synthesis",
"kBasalSynPIP2_PIP2_synthesis",
"PIP_basal_PIP_synthesis",
"kBasalSynPIP_PIP_synthesis",
"kStimSynPIP_PIP_synthesis",
"tauPIPsyn_PIP_synthesis",
"PIPsyndecay_PIP_synthesis",
"Kon_reaction2",
"dinh_reaction2",
"dact_Ca_binds_IP3R",
"PIP2_synthesis_I",
"Shc_phosphorylation_Km",
"Shc_phosphorylation_I",
"Ca_bind_CAPRI_Kr",
"Ca_bind_CAPRI_Kf",
"IP3_degradation_IP3_basal",
"IP3_degradation_kIP3deg",
"RasGTP_depal_translocate_Kr",
"RasGTP_depal_translocate_Kf",
"RasGTP_depal_translocate_I",
"ca_bind_rasGRP_Kr",
"ca_bind_rasGRP_Kf",
"RasGRP_DAG_GEF_Km",
"RasGRP_DAG_GEF_I",
"rasGTP_pal_translocation_Kr",
"rasGTP_pal_translocation_Kf",
"rasGTP_pal_translocation_I",
"PLCg_dephos_Kr",
"PLCg_dephos_Kf",
"PLCg_dephos_I",
"basal_GAP_Kr",
"basal_GAP_Kf",
"basal_GAP_I",
"CAPRI_translocation_Kr",
"CAPRI_translocation_Kf",
"CAPRI_translocation_I",
"reaction5_Kr",
"reaction5_Kf",
"reaction5_I",
"RasGDP_pal_Kr",
"RasGDP_pal_Kf",
"RasGDP_pal_I",
"CAPRI_GAP_Km",
"CAPRI_GAP_I",
"RasGDPpal_translocation_Kr",
"RasGDPpal_translocation_Kf",
"RasGDPpal_translocation_I",
"sos_grb2_binding_Kr",
"sos_grb2_binding_Kf",
"RasGDP_depal_translocate_Kr",
"RasGDP_depal_translocate_Kf",
"RasGDP_depal_translocate_I",
"Ras_GTP_palm1_Kr",
"Ras_GTP_palm1_Kf",
"Ras_GTP_palm1_I",
"RasPal_basal_GAP_Kr",
"RasPal_basal_GAP_Kf",
"basal_cyt_depal_GEF_Kr",
"basal_cyt_depal_GEF_Kf",
"caPLCe_gen_DAG_kact",
"caPLCe_gen_DAG_I",
"CaRasGRP_act_RasGM_Km",
"CaRasGRP_act_RasGM_I",
"PIP2_hydrolysis_k_PIP2hyd",
"PIP2_hydrolysis_I",
"Sos_activation_Kr",
"Sos_activation_Kf",
"Sos_activation_I",
"PIP_synthesis_I",
"EGF_internalization_Kr",
"EGF_internalization_Kf",
"EGF_internalization_I",
"calcium_buffer_Kr",
"calcium_buffer_Kf",
"RasGM_basal_GAP_Vmax",
"RasGM_basal_GAP_Km",
"RasGM_basal_GAP_I",
"reaction0_Kr",
"reaction0_Kf",
"reaction0_I",
"EGFR_binding_Kr",
"EGFR_binding_Kf",
"EGFR_binding_I",
"ca_act_PLCe_Kr",
"ca_act_PLCe_Kf",
"ras_act_PLCe_Kr",
"ras_act_PLCe_Kf",
"ras_act_PLCe_I",
"RasGDP_depal2_Kr",
"RasGDP_depal2_Kf",
"RasGDP_depal2_I",
"CaRasGRP_translocation_Kr",
"CaRasGRP_translocation_Kf",
"CaRasGRP_translocation_I",
"reaction2_I",
"EGF_act_PLCgamma_Km",
"EGF_act_PLCgamma_I",
"Ca_binds_IP3R_Kf",
"Ca_binds_IP3R_I",
"reaction7_Kr",
"reaction7_Kf",
"reaction7_I",
"Sos_act_RasPM_Km",
"Sos_act_RasPM_I",
"flux1_kP",
"flux1_vP",
"flux1_I",
"flux0_singleChanFlux",
"flux0_dI",
"flux0_I",
"flux2_vL",
"flux2_I",
"Ras_PM_depal1_Kr",
"Ras_PM_depal1_Kf",
"Ras_PM_depal1_I",};
    }

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    virtual std::vector<std::string> getStateIds() const override {
        return std::vector<std::string>{"RasGTP_Golgi_GM",
"EGF_EC",
"CAPRI_cyt",
"serca",
"PIP_PM",
"PI_PM",
"Shc_PM",
"CaCAPRI_PM_PM",
"RactCa",
"Shc_star_PM",
"EGFR_PM",
"PLC_act_PM",
"RasGTP_pal_cyt",
"PLC_PM",
"PIP2_PM",
"Activated_EGFR_PM",
"Ca",
"Ract",
"Rinh",
"RinhCa",
"IP3",
"RasGDP_Golgi_GM",
"Ca_RasGRP_GM_GM",
"DAG_GM_GM",
"RasGRP_DAG_GM",
"CaCAPRI_cyt",
"DAG_PM",
"RasGTP_depal_cyt",
"RasGDP_depal_cyt",
"RasGDP_pal_cyt",
"Ca_PLCe_cyt",
"Ras_CaPLCe_GM",
"PIP2_GM_GM",
"ER_erMembrane",
"Ca_ER",
"Sos_cyt",
"Grb2_cyt",
"PLCe_cyt",
"buffer_cyt",
"ca_buffer_cyt",
"SosGrb2_cyt",
"SGS_PM",
"RasGTP_PM",
"RasGDP_PM",
"RasGRP_cyt",
"CaRasGRP1_cyt",};
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
"y45",};
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
