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
extern void J_Levchenko2000a(realtype *J, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w,
                            const realtype *dwdx);
extern void JB_Levchenko2000a(realtype *JB, const realtype t, const realtype *x,
                             const realtype *p, const realtype *k,
                             const realtype *h, const realtype *xB,
                             const realtype *w, const realtype *dwdx);
extern void JDiag_Levchenko2000a(realtype *JDiag, const realtype t,
                                const realtype *x, const realtype *p,
                                const realtype *k, const realtype *h,
                                const realtype *w, const realtype *dwdx);
extern void JSparse_Levchenko2000a(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JSparse_colptrs_Levchenko2000a(sunindextype *colptrs);
extern void JSparse_rowvals_Levchenko2000a(sunindextype *rowvals);
extern void JSparseB_Levchenko2000a(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void JSparseB_colptrs_Levchenko2000a(sunindextype *colptrs);
extern void JSparseB_rowvals_Levchenko2000a(sunindextype *rowvals);
extern void Jy_Levchenko2000a(realtype *nllh, const int iy, const realtype *p,
                             const realtype *k, const realtype *y,
                             const realtype *sigmay, const realtype *my);
extern void dJydsigmay_Levchenko2000a(realtype *dJydsigmay, const int iy,
                                     const realtype *p, const realtype *k,
                                     const realtype *y, const realtype *sigmay,
                                     const realtype *my);
extern void dJydy_Levchenko2000a(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_Levchenko2000a(sunindextype *colptrs, int index);
extern void dJydy_rowvals_Levchenko2000a(sunindextype *rowvals, int index);
extern void dwdp_Levchenko2000a(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip);
extern void dwdx_Levchenko2000a(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl);
extern void dwdx_colptrs_Levchenko2000a(sunindextype *colptrs);
extern void dwdx_rowvals_Levchenko2000a(sunindextype *rowvals);
extern void dxdotdw_Levchenko2000a(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdw_colptrs_Levchenko2000a(sunindextype *colptrs);
extern void dxdotdw_rowvals_Levchenko2000a(sunindextype *rowvals);
extern void dxdotdp_Levchenko2000a(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w);
extern void dydx_Levchenko2000a(realtype *dydx, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w, const realtype *dwdx);
extern void dydp_Levchenko2000a(realtype *dydp, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const int ip, const realtype *w,
                               const realtype *dwp);
extern void dsigmaydp_Levchenko2000a(realtype *dsigmaydp, const realtype t,
                                    const realtype *p, const realtype *k,
                                    const int ip);
extern void sigmay_Levchenko2000a(realtype *sigmay, const realtype t,
                                 const realtype *p, const realtype *k);
extern void w_Levchenko2000a(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl);
extern void x0_Levchenko2000a(realtype *x0, const realtype t, const realtype *p,
                             const realtype *k);
extern void x0_fixedParameters_Levchenko2000a(realtype *x0, const realtype t,
                                             const realtype *p,
                                             const realtype *k);
extern void sx0_Levchenko2000a(realtype *sx0, const realtype t,
                              const realtype *x0, const realtype *p,
                              const realtype *k, const int ip);
extern void sx0_fixedParameters_Levchenko2000a(realtype *sx0, const realtype t,
                                              const realtype *x0,
                                              const realtype *p,
                                              const realtype *k, const int ip);
extern void xdot_Levchenko2000a(realtype *xdot, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w);
extern void y_Levchenko2000a(realtype *y, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w);

extern void x_solver_Levchenko2000a(realtype *x_solver, const realtype *x_rdata);
extern void total_cl_Levchenko2000a(realtype *total_cl, const realtype *x_rdata);

/**
 * @brief AMICI-generated model subclass.
 */
class Model_Levchenko2000a : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_Levchenko2000a()
        : amici::Model_ODE(
              86,                                // nx_rdata
              86,                            // nxtrue_rdata
              86,                               // nx_solver
              86,                           // nxtrue_solver
              86,                                      // ny
              86,                                  // nytrue
              0,                                      // nz
              0,                                  // nztrue
              0,                                  // nevent
              1,                              // nobjective
              300,                                      // nw
              430,                                   // ndwdx
              90000,                                   // ndwdp
              886,                                // ndxdotdw
              std::vector<int>{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},                                  // ndjydy
              906,                                     // nnz
              86,                                     // ubw
              86,                                     // lbw
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{1.0, 0.4, 0.1, 0.5, 0.5, 0.1, 3.3, 0.42, 0.1, 10.0, 0.8, 0.1, 3.3, 0.4, 0.1, 10.0, 0.8, 0.1, 20.0, 0.6, 0.1, 5.0, 0.4, 0.1, 20.0, 0.6, 0.1, 5.0, 0.4, 0.1, 10.0, 0.5, 10.0, 0.5, 10.0, 0.5, 10.0, 0.5, 10.0, 0.5, 10.0, 0.5, 10.0, 0.5, 10.0, 0.5, 10.0, 0.5, 10.0, 0.5, 10.0, 0.5, 10.0, 0.5, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 10.0, 0.5, 10.0, 0.5, 10.0, 0.5, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 10.0, 0.5, 10.0, 0.5, 10.0, 0.5, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 10.0, 0.5, 10.0, 0.5, 10.0, 0.5, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 10.0, 0.5, 10.0, 0.5, 10.0, 0.5, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 0.0, 0.05, 10.0, 0.5, 0.0, 0.05, 10.0, 0.5, 0.0, 0.05, 10.0, 0.5, 0.0, 0.05, 10.0, 0.5, 0.0, 0.05, 10.0, 0.5, 0.0, 0.05, 10.0, 0.5, 0.0, 0.05, 10.0, 0.5, 0.0, 0.05, 10.0, 0.5, 0.0, 0.05, 10.0, 0.5, 0.0, 0.05, 10.0, 0.5, 0.0, 0.05, 10.0, 0.5, 0.0, 0.05, 10.0, 0.5, 0.0, 0.05, 10.0, 0.5, 0.0, 0.05, 10.0, 0.5, 0.0, 0.05, 10.0, 0.5, 0.0, 0.05, 10.0, 0.5, 0.0, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 100.0, 0.0, 0.1, 100.0, 0.0, 0.1, 100.0, 0.0, 0.1, 100.0, 0.0, 0.1, 100.0, 0.0, 0.1, 100.0, 0.0, 0.1, 100.0, 0.0, 0.1, 100.0, 0.0, 0.1, 100.0, 0.0, 0.1, 100.0, 0.0, 0.1, 100.0, 0.0, 0.1, 100.0, 0.0, 0.1, 100.0, 0.0, 0.1, 100.0, 0.0, 0.1, 100.0, 0.0, 0.1, 100.0, 0.0, 0.1},       // dynamic parameters
              std::vector<realtype>{}, // fixedParameters
              std::vector<int>{},                          // plist
              std::vector<realtype>(86, 0.0),   // idlist
              std::vector<int>{}                           // z2event
          ) {}

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    virtual amici::Model *clone() const override {
        return new Model_Levchenko2000a(*this);
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
        J_Levchenko2000a(J, t, x, p, k, h, w, dwdx);
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
        JB_Levchenko2000a(JB, t, x, p, k, h, xB, w, dwdx);
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
        JDiag_Levchenko2000a(JDiag, t, x, p, k, h, w, dwdx);
    }

    virtual void fJSparse(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        JSparse_Levchenko2000a( JSparse,  t,  x,  p,  k,  h,  w,  dwdx);
    }


    virtual void fJSparse_colptrs(sunindextype *colptrs) override {
        JSparse_colptrs_Levchenko2000a(colptrs);
    }


    virtual void fJSparse_rowvals(sunindextype *rowvals) override {
        JSparse_rowvals_Levchenko2000a(rowvals);
    }


    virtual void fJSparseB(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        JSparseB_Levchenko2000a( JSparseB,  t,  x,  p,  k,  h,  xB,  w,  dwdx);
    }


    virtual void fJSparseB_colptrs(sunindextype *colptrs) override {
        JSparseB_colptrs_Levchenko2000a(colptrs);
    }


    virtual void fJSparseB_rowvals(sunindextype *rowvals) override {
        JSparseB_rowvals_Levchenko2000a(rowvals);
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
        Jy_Levchenko2000a(nllh, iy, p, k, y, sigmay, my);
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
        dJydsigmay_Levchenko2000a(dJydsigma, iy, p, k, y, sigmay, my);
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
        dsigmaydp_Levchenko2000a(dsigmaydp, t, p, k, ip);
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
        dJydy_Levchenko2000a( dJydy,  iy,  p,  k,  y,  sigmay,  my);
    }

    virtual void fdJydy_colptrs(sunindextype *colptrs, int index) override {
        dJydy_colptrs_Levchenko2000a(colptrs, index);
    }

    virtual void fdJydy_rowvals(sunindextype *rowvals, int index) override {
        dJydy_rowvals_Levchenko2000a(rowvals, index);
    }


    virtual void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip) override {
        dwdp_Levchenko2000a( dwdp,  t,  x,  p,  k,  h,  w,  tcl,  dtcldp,  ip);
    }


    virtual void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {
        dwdx_Levchenko2000a( dwdx,  t,  x,  p,  k,  h,  w,  tcl);
    }


    virtual void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdw_Levchenko2000a( dxdotdw,  t,  x,  p,  k,  h,  w);
    }


    virtual void fdxdotdw_colptrs(sunindextype *colptrs) override {
        dxdotdw_colptrs_Levchenko2000a(colptrs);
    }


    virtual void fdxdotdw_rowvals(sunindextype *rowvals) override {
        dxdotdw_rowvals_Levchenko2000a(rowvals);
    }


    virtual void fdxdotdp(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w) override {
        dxdotdp_Levchenko2000a( dxdotdp,  t,  x,  p,  k,  h,  ip,  w);
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
        dydx_Levchenko2000a(dydx, t, x, p, k, h, w, dwdx);
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
        dydp_Levchenko2000a(dydp, t, x, p, k, h, ip, w, dwdp);
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
        sigmay_Levchenko2000a(sigmay, t, p, k);
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
        sx0_Levchenko2000a(sx0, t, x0, p, k, ip);
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
        sx0_fixedParameters_Levchenko2000a(sx0, t, x0, p, k, ip);
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
        w_Levchenko2000a( w,  t,  x,  p,  k,  h,  tcl);
    }


    /** model specific implementation of fx0
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0(realtype *x0, const realtype t, const realtype *p,
                     const realtype *k) override {
        x0_Levchenko2000a(x0, t, p, k);
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
        x0_fixedParameters_Levchenko2000a(x0, t, p, k);
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
        xdot_Levchenko2000a(xdot, t, x, p, k, h, w);
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
        y_Levchenko2000a(y, t, x, p, k, h, w);
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
        x_solver_Levchenko2000a( x_solver,  x_rdata);
    }


    virtual void ftotal_cl(realtype *total_cl, const realtype *x_rdata) override {
        total_cl_Levchenko2000a( total_cl,  x_rdata);
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
        return std::vector<std::string>{"MAPK phosphatase",
"MEK phosphatase",
"RAF kinase",
"RAF phosphatase",
"MAPK",
"MAPK-P",
"MAPK-PP",
"MEK",
"MEK-P",
"MEK-PP",
"RAF",
"RAF-P",
"MAPK_MEK-PP",
"MAPK-P_MEK-PP",
"MEK_RAF-P",
"MEK-P_RAF-P",
"MAPK-P_MAPKPase",
"MAPK-PP_MAPKPase",
"MEK-P_MEKPase",
"MEK-PP_MEKPase",
"RAF_RAFK",
"RAF-P_RAFPase",
"Scaffold",
"Scaffold_RAF",
"Scaffold_RAF-P",
"Scaffold_MEK",
"Scaffold_MEK_RAF",
"Scaffold_MEK_RAF-P",
"Scaffold_MEK-P",
"Scaffold_MEK-P_RAF",
"Scaffold_MEK-P_RAF-P",
"Scaffold_MEK-PP",
"Scaffold_MEK-PP_RAF",
"Scaffold_MEK-PP_RAF-P",
"Scaffold_MAPK",
"Scaffold_MAPK_RAF",
"Scaffold_MAPK_RAF-P",
"Scaffold_MAPK_MEK",
"Scaffold_MAPK_MEK_RAF",
"Scaffold_MAPK_MEK_RAF-P",
"Scaffold_MAPK_MEK-P",
"Scaffold_MAPK_MEK-P_RAF",
"Scaffold_MAPK_MEK-P_RAF-P",
"Scaffold_MAPK_MEK-PP",
"Scaffold_MAPK_MEK-PP_RAF",
"Scaffold_MAPK_MEK-PP_RAF-P",
"Scaffold_MAPK-P",
"Scaffold_MAPK-P_RAF",
"Scaffold_MAPK-P_RAF-P",
"Scaffold_MAPK-P_MEK",
"Scaffold_MAPK-P_MEK_RAF",
"Scaffold_MAPK-P_MEK_RAF-P",
"Scaffold_MAPK-P_MEK-P",
"Scaffold_MAPK-P_MEK-P_RAF",
"Scaffold_MAPK-P_MEK-P_RAF-P",
"Scaffold_MAPK-P_MEK-PP",
"Scaffold_MAPK-P_MEK-PP_RAF",
"Scaffold_MAPK-P_MEK-PP_RAF-P",
"Scaffold_MAPK-PP",
"Scaffold_MAPK-PP_RAF",
"Scaffold_MAPK-PP_RAF-P",
"Scaffold_MAPK-PP_MEK",
"Scaffold_MAPK-PP_MEK_RAF",
"Scaffold_MAPK-PP_MEK_RAF-P",
"Scaffold_MAPK-PP_MEK-P",
"Scaffold_MAPK-PP_MEK-P_RAF",
"Scaffold_MAPK-PP_MEK-P_RAF-P",
"Scaffold_MAPK-PP_MEK-PP",
"Scaffold_MAPK-PP_MEK-PP_RAF",
"Scaffold_MAPK-PP_MEK-PP_RAF-P",
"Scaffold_RAF",
"Scaffold_MEK_RAF",
"Scaffold_MEK-P_RAF",
"Scaffold_MEK-PP_RAF",
"Scaffold_MAPK_RAF",
"Scaffold_MAPK_MEK_RAF",
"Scaffold_MAPK_MEK-P_RAF",
"Scaffold_MAPK_MEK-PP_RAF",
"Scaffold_MAPK-P_RAF",
"Scaffold_MAPK-P_MEK_RAF",
"Scaffold_MAPK-P_MEK-P_RAF",
"Scaffold_MAPK-P_MEK-PP_RAF",
"Scaffold_MAPK-PP_RAF",
"Scaffold_MAPK-PP_MEK_RAF",
"Scaffold_MAPK-PP_MEK-P_RAF",
"Scaffold_MAPK-PP_MEK-PP_RAF",};
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
"x85",};
    }

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getParameterIds() const override {
        return std::vector<std::string>{"Reaction1_a1",
"Reaction2_d1",
"Reaction3_k1",
"Reaction4_a2",
"Reaction5_d2",
"Reaction6_k2",
"Reaction7_a3",
"Reaction8_d3",
"Reaction9_k3",
"Reaction10_a4",
"Reaction11_d4",
"Reaction12_k4",
"Reaction13_a5",
"Reaction14_d5",
"Reaction15_k5",
"Reaction16_a6",
"Reaction17_d6",
"Reaction18_k6",
"Reaction19_a7",
"Reaction20_d7",
"Reaction21_k7",
"Reaction22_a8",
"Reaction23_d8",
"Reaction24_k8",
"Reaction25_a9",
"Reaction26_d9",
"Reaction27_k9",
"Reaction28_a10",
"Reaction29_d10",
"Reaction30_k10",
"Reaction31_kon",
"Reaction32_koff",
"Reaction33_kon",
"Reaction34_koff",
"Reaction35_kon",
"Reaction36_koff",
"Reaction37_kon",
"Reaction38_koff",
"Reaction39_kon",
"Reaction40_koff",
"Reaction41_kon",
"Reaction42_koff",
"Reaction43_kon",
"Reaction44_koff",
"Reaction45_kon",
"Reaction46_koff",
"Reaction47_kon",
"Reaction48_koff",
"Reaction49_kon",
"Reaction50_koff",
"Reaction51_kon",
"Reaction52_koff",
"Reaction53_kon",
"Reaction54_koff",
"Reaction55_kpon",
"Reaction56_kpoff",
"Reaction57_kpon",
"Reaction58_kpoff",
"Reaction59_kpon",
"Reaction60_kpoff",
"Reaction61_kpon",
"Reaction62_kpoff",
"Reaction63_kpon",
"Reaction64_kpoff",
"Reaction65_kpon",
"Reaction66_kpoff",
"Reaction67_kpon",
"Reaction68_kpoff",
"Reaction69_kpon",
"Reaction70_kpoff",
"Reaction71_kpon",
"Reaction72_kpoff",
"Reaction73_kpon",
"Reaction74_kpoff",
"Reaction75_kpon",
"Reaction76_kpoff",
"Reaction77_kpon",
"Reaction78_kpoff",
"Reaction79_kpon",
"Reaction80_kpoff",
"Reaction81_kpon",
"Reaction82_kpoff",
"Reaction83_kpon",
"Reaction84_kpoff",
"Reaction85_kpon",
"Reaction86_kpoff",
"Reaction87_kpon",
"Reaction88_kpoff",
"Reaction89_kpon",
"Reaction90_kpoff",
"Reaction91_kpon",
"Reaction92_kpoff",
"Reaction93_kpon",
"Reaction94_kpoff",
"Reaction95_kpon",
"Reaction96_kpoff",
"Reaction97_kpon",
"Reaction98_kpoff",
"Reaction99_kpon",
"Reaction100_kpoff",
"Reaction101_kpon",
"Reaction102_kpoff",
"Reaction103_kon",
"Reaction104_koff",
"Reaction105_kon",
"Reaction106_koff",
"Reaction107_kon",
"Reaction108_koff",
"Reaction109_kpon",
"Reaction110_kpoff",
"Reaction111_kpon",
"Reaction112_kpoff",
"Reaction113_kpon",
"Reaction114_kpoff",
"Reaction115_kpon",
"Reaction116_kpoff",
"Reaction117_kpon",
"Reaction118_kpoff",
"Reaction119_kpon",
"Reaction120_kpoff",
"Reaction121_kon",
"Reaction122_koff",
"Reaction123_kon",
"Reaction124_koff",
"Reaction125_kon",
"Reaction126_koff",
"Reaction127_kpon",
"Reaction128_kpoff",
"Reaction129_kpon",
"Reaction130_kpoff",
"Reaction131_kpon",
"Reaction132_kpoff",
"Reaction133_kpon",
"Reaction134_kpoff",
"Reaction135_kpon",
"Reaction136_kpoff",
"Reaction137_kpon",
"Reaction138_kpoff",
"Reaction139_kon",
"Reaction140_koff",
"Reaction141_kon",
"Reaction142_koff",
"Reaction143_kon",
"Reaction144_koff",
"Reaction145_kpon",
"Reaction146_kpoff",
"Reaction147_kpon",
"Reaction148_kpoff",
"Reaction149_kpon",
"Reaction150_kpoff",
"Reaction151_kpon",
"Reaction152_kpoff",
"Reaction153_kpon",
"Reaction154_kpoff",
"Reaction155_kpon",
"Reaction156_kpoff",
"Reaction157_kon",
"Reaction158_koff",
"Reaction159_kon",
"Reaction160_koff",
"Reaction161_kon",
"Reaction162_koff",
"Reaction163_kpon",
"Reaction164_kpoff",
"Reaction165_kpon",
"Reaction166_kpoff",
"Reaction167_kpon",
"Reaction168_kpoff",
"Reaction169_kpon",
"Reaction170_kpoff",
"Reaction171_kpon",
"Reaction172_kpoff",
"Reaction173_kpon",
"Reaction174_kpoff",
"Reaction175_kon",
"Reaction176_koff",
"Reaction177_kpon",
"Reaction178_kpoff",
"Reaction179_kon",
"Reaction180_koff",
"Reaction181_kpon",
"Reaction182_kpoff",
"Reaction183_kon",
"Reaction184_koff",
"Reaction185_kpon",
"Reaction186_kpoff",
"Reaction187_kon",
"Reaction188_koff",
"Reaction189_kpon",
"Reaction190_kpoff",
"Reaction191_kon",
"Reaction192_koff",
"Reaction193_kpon",
"Reaction194_kpoff",
"Reaction195_kon",
"Reaction196_koff",
"Reaction197_kpon",
"Reaction198_kpoff",
"Reaction199_kon",
"Reaction200_koff",
"Reaction201_kpon",
"Reaction202_kpoff",
"Reaction203_kon",
"Reaction204_koff",
"Reaction205_kpon",
"Reaction206_kpoff",
"Reaction207_kon",
"Reaction208_koff",
"Reaction209_kpon",
"Reaction210_kpoff",
"Reaction211_kon",
"Reaction212_koff",
"Reaction213_kpon",
"Reaction214_kpoff",
"Reaction215_kon",
"Reaction216_koff",
"Reaction217_kpon",
"Reaction218_kpoff",
"Reaction219_kon",
"Reaction220_koff",
"Reaction221_kpon",
"Reaction222_kpoff",
"Reaction223_kon",
"Reaction224_koff",
"Reaction225_kpon",
"Reaction226_kpoff",
"Reaction227_kon",
"Reaction228_koff",
"Reaction229_kpon",
"Reaction230_kpoff",
"Reaction231_kon",
"Reaction232_koff",
"Reaction233_kpon",
"Reaction234_kpoff",
"Reaction235_kon",
"Reaction236_koff",
"Reaction237_kpon",
"Reaction238_kpoff",
"Reaction239_k7",
"Reaction240_k7",
"Reaction241_k7",
"Reaction242_k9a",
"Reaction243_k9a",
"Reaction244_k9a",
"Reaction245_k3",
"Reaction246_k5a",
"Reaction247_k3",
"Reaction248_k5a",
"Reaction249_k3",
"Reaction250_k5a",
"Reaction251_k3",
"Reaction252_k5a",
"Reaction253_k1a",
"Reaction254_d1a",
"Reaction255_k1",
"Reaction256_k1a",
"Reaction257_d1a",
"Reaction258_k1",
"Reaction259_k1a",
"Reaction260_d1a",
"Reaction261_k1",
"Reaction262_k1a",
"Reaction263_d1a",
"Reaction264_k1",
"Reaction265_k1a",
"Reaction266_d1a",
"Reaction267_k1",
"Reaction268_k1a",
"Reaction269_d1a",
"Reaction270_k1",
"Reaction271_k1a",
"Reaction272_d1a",
"Reaction273_k1",
"Reaction274_k1a",
"Reaction275_d1a",
"Reaction276_k1",
"Reaction277_k1a",
"Reaction278_d1a",
"Reaction279_k1",
"Reaction280_k1a",
"Reaction281_d1a",
"Reaction282_k1",
"Reaction283_k1a",
"Reaction284_d1a",
"Reaction285_k1",
"Reaction286_k1a",
"Reaction287_d1a",
"Reaction288_k1",
"Reaction289_k1a",
"Reaction290_d1a",
"Reaction291_k1",
"Reaction292_k1a",
"Reaction293_d1a",
"Reaction294_k1",
"Reaction295_k1a",
"Reaction296_d1a",
"Reaction297_k1",
"Reaction298_k1a",
"Reaction299_d1a",
"Reaction300_k1",};
    }

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    virtual std::vector<std::string> getStateIds() const override {
        return std::vector<std::string>{"MAPKP",
"MEKP",
"RAFK",
"RAFP",
"K_1_0",
"K_1_1",
"K_1_2",
"K_2_0",
"K_2_1",
"K_2_2",
"K_3_0",
"K_3_1",
"K_K_1_0_2_2",
"K_K_1_1_2_2",
"K_K_2_0_3_1",
"K_K_2_1_3_1",
"K_MAPKP_1_1",
"K_MAPKP_1_2",
"K_MEKP_2_1",
"K_MEKP_2_2",
"K_RAFK_3_0",
"K_RAFP_3_1",
"S_m1_m1_m1",
"S_m1_m1_0",
"S_m1_m1_1",
"S_m1_0_m1",
"S_m1_0_0",
"S_m1_0_1",
"S_m1_1_m1",
"S_m1_1_0",
"S_m1_1_1",
"S_m1_2_m1",
"S_m1_2_0",
"S_m1_2_1",
"S_0_m1_m1",
"S_0_m1_0",
"S_0_m1_1",
"S_0_0_m1",
"S_0_0_0",
"S_0_0_1",
"S_0_1_m1",
"S_0_1_0",
"S_0_1_1",
"S_0_2_m1",
"S_0_2_0",
"S_0_2_1",
"S_1_m1_m1",
"S_1_m1_0",
"S_1_m1_1",
"S_1_0_m1",
"S_1_0_0",
"S_1_0_1",
"S_1_1_m1",
"S_1_1_0",
"S_1_1_1",
"S_1_2_m1",
"S_1_2_0",
"S_1_2_1",
"S_2_m1_m1",
"S_2_m1_0",
"S_2_m1_1",
"S_2_0_m1",
"S_2_0_0",
"S_2_0_1",
"S_2_1_m1",
"S_2_1_0",
"S_2_1_1",
"S_2_2_m1",
"S_2_2_0",
"S_2_2_1",
"S_RAFK_m1_m1_0",
"S_RAFK_m1_0_0",
"S_RAFK_m1_1_0",
"S_RAFK_m1_2_0",
"S_RAFK_0_m1_0",
"S_RAFK_0_0_0",
"S_RAFK_0_1_0",
"S_RAFK_0_2_0",
"S_RAFK_1_m1_0",
"S_RAFK_1_0_0",
"S_RAFK_1_1_0",
"S_RAFK_1_2_0",
"S_RAFK_2_m1_0",
"S_RAFK_2_0_0",
"S_RAFK_2_1_0",
"S_RAFK_2_2_0",};
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
"y85",};
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
