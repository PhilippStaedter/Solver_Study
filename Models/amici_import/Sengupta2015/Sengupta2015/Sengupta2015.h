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
extern void J_Sengupta2015(realtype *J, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w,
                            const realtype *dwdx);
extern void JB_Sengupta2015(realtype *JB, const realtype t, const realtype *x,
                             const realtype *p, const realtype *k,
                             const realtype *h, const realtype *xB,
                             const realtype *w, const realtype *dwdx);
extern void JDiag_Sengupta2015(realtype *JDiag, const realtype t,
                                const realtype *x, const realtype *p,
                                const realtype *k, const realtype *h,
                                const realtype *w, const realtype *dwdx);
extern void JSparse_Sengupta2015(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JSparse_colptrs_Sengupta2015(sunindextype *colptrs);
extern void JSparse_rowvals_Sengupta2015(sunindextype *rowvals);
extern void JSparseB_Sengupta2015(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void JSparseB_colptrs_Sengupta2015(sunindextype *colptrs);
extern void JSparseB_rowvals_Sengupta2015(sunindextype *rowvals);
extern void Jy_Sengupta2015(realtype *nllh, const int iy, const realtype *p,
                             const realtype *k, const realtype *y,
                             const realtype *sigmay, const realtype *my);
extern void dJydsigmay_Sengupta2015(realtype *dJydsigmay, const int iy,
                                     const realtype *p, const realtype *k,
                                     const realtype *y, const realtype *sigmay,
                                     const realtype *my);
extern void dJydy_Sengupta2015(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_Sengupta2015(sunindextype *colptrs, int index);
extern void dJydy_rowvals_Sengupta2015(sunindextype *rowvals, int index);
extern void dwdp_Sengupta2015(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip);
extern void dwdx_Sengupta2015(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl);
extern void dwdx_colptrs_Sengupta2015(sunindextype *colptrs);
extern void dwdx_rowvals_Sengupta2015(sunindextype *rowvals);
extern void dxdotdw_Sengupta2015(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdw_colptrs_Sengupta2015(sunindextype *colptrs);
extern void dxdotdw_rowvals_Sengupta2015(sunindextype *rowvals);
extern void dxdotdp_Sengupta2015(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w);
extern void dydx_Sengupta2015(realtype *dydx, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w, const realtype *dwdx);
extern void dydp_Sengupta2015(realtype *dydp, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const int ip, const realtype *w,
                               const realtype *dwp);
extern void dsigmaydp_Sengupta2015(realtype *dsigmaydp, const realtype t,
                                    const realtype *p, const realtype *k,
                                    const int ip);
extern void sigmay_Sengupta2015(realtype *sigmay, const realtype t,
                                 const realtype *p, const realtype *k);
extern void w_Sengupta2015(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl);
extern void x0_Sengupta2015(realtype *x0, const realtype t, const realtype *p,
                             const realtype *k);
extern void x0_fixedParameters_Sengupta2015(realtype *x0, const realtype t,
                                             const realtype *p,
                                             const realtype *k);
extern void sx0_Sengupta2015(realtype *sx0, const realtype t,
                              const realtype *x0, const realtype *p,
                              const realtype *k, const int ip);
extern void sx0_fixedParameters_Sengupta2015(realtype *sx0, const realtype t,
                                              const realtype *x0,
                                              const realtype *p,
                                              const realtype *k, const int ip);
extern void xdot_Sengupta2015(realtype *xdot, const realtype t,
                               const realtype *x, const realtype *p,
                               const realtype *k, const realtype *h,
                               const realtype *w);
extern void y_Sengupta2015(realtype *y, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k,
                            const realtype *h, const realtype *w);

extern void x_solver_Sengupta2015(realtype *x_solver, const realtype *x_rdata);
extern void total_cl_Sengupta2015(realtype *total_cl, const realtype *x_rdata);

/**
 * @brief AMICI-generated model subclass.
 */
class Model_Sengupta2015 : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_Sengupta2015()
        : amici::Model_ODE(
              240,                                // nx_rdata
              240,                            // nxtrue_rdata
              240,                               // nx_solver
              240,                           // nxtrue_solver
              240,                                      // ny
              240,                                  // nytrue
              0,                                      // nz
              0,                                  // nztrue
              0,                                  // nevent
              1,                              // nobjective
              175,                                      // nw
              195,                                   // ndwdx
              47600,                                   // ndwdp
              545,                                // ndxdotdw
              std::vector<int>{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},                                  // ndjydy
              571,                                     // nnz
              240,                                     // ubw
              240,                                     // lbw
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.069, 1.0, 0.048, 1.0, 0.04, 1.0, 1.0, 875.0, 1.0, 0.013, 1.0, 0.00038, 1.0, 1.3, 1.0, 13.6, 1.0, 13.0, 1.0, 0.1, 1.0, 5.8, 1.0, 0.11, 1.0, 73.0, 1.0, 490.0, 1.0, 0.09, 1.0, 2.0, 1.0, 0.000213, 1.0, 33.0, 1.0, 300.0, 1.0, 300.0, 1.0, 791.0, 1.0, 0.11, 1.0, 40.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.04, 1.0, 18.2, 1.0, 0.157, 1.0, 13.6, 1.0, 140.0, 1.0, 0.066, 1.0, 0.13, 1.0, 0.045, 1.0, 0.1, 1.0, 0.08, 1.0, 0.1, 1.0, 10.0, 1.0, 0.58, 1.0, 1.7, 1.0, 0.31, 1.0, 1.37, 1.0, 970.0, 1.0, 0.13, 1.0, 294.0, 1.0, 1.28, 1.0, 12.6, 1.0, 18.7, 1.0, 0.08, 1.0, 0.19, 1.0, 0.19, 1.0, 2900.0, 1.0, 2900.0, 1.0, 16.0, 1.0, 16.0, 1.0, 16.0, 1.0, 16.0, 1.0, 34.5, 1.0, 16.0, 1.0, 16.0, 1.0, 16.0, 1.0, 16.0, 1.0, 16.0, 1.0, 16.0, 1.0, 100.0, 1.0, 16.0, 1.0, 3.0, 1.0, 3.0, 1.0, 3.0, 1.0, 3.0, 1.0, 3.0, 1.0, 3.0, 1.0, 3.0, 1.0, 3.0, 1.0, 3.0, 1.0, 3.0, 1.0, 100.0, 1.0, 100.0, 1.0, 100.0, 1.0, 100.0, 1.0, 100.0, 1.0, 100.0, 1.0, 100.0, 1.0, 100.0, 1.0, 100.0, 1.0, 100.0, 1.0, 34.5, 1.0, 34.5, 1.0, 34.5, 1.0, 34.5, 1.0, 34.5, 1.0, 34.5, 1.0, 34.5, 1.0, 34.5, 1.0, 34.5, 1.0, 34.5, 1.0, 1.16, 1.0, 9.6, 1.0},       // dynamic parameters
              std::vector<realtype>{}, // fixedParameters
              std::vector<int>{},                          // plist
              std::vector<realtype>(240, 0.0),   // idlist
              std::vector<int>{}                           // z2event
          ) {}

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    virtual amici::Model *clone() const override {
        return new Model_Sengupta2015(*this);
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
        J_Sengupta2015(J, t, x, p, k, h, w, dwdx);
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
        JB_Sengupta2015(JB, t, x, p, k, h, xB, w, dwdx);
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
        JDiag_Sengupta2015(JDiag, t, x, p, k, h, w, dwdx);
    }

    virtual void fJSparse(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        JSparse_Sengupta2015( JSparse,  t,  x,  p,  k,  h,  w,  dwdx);
    }


    virtual void fJSparse_colptrs(sunindextype *colptrs) override {
        JSparse_colptrs_Sengupta2015(colptrs);
    }


    virtual void fJSparse_rowvals(sunindextype *rowvals) override {
        JSparse_rowvals_Sengupta2015(rowvals);
    }


    virtual void fJSparseB(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        JSparseB_Sengupta2015( JSparseB,  t,  x,  p,  k,  h,  xB,  w,  dwdx);
    }


    virtual void fJSparseB_colptrs(sunindextype *colptrs) override {
        JSparseB_colptrs_Sengupta2015(colptrs);
    }


    virtual void fJSparseB_rowvals(sunindextype *rowvals) override {
        JSparseB_rowvals_Sengupta2015(rowvals);
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
        Jy_Sengupta2015(nllh, iy, p, k, y, sigmay, my);
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
        dJydsigmay_Sengupta2015(dJydsigma, iy, p, k, y, sigmay, my);
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
        dsigmaydp_Sengupta2015(dsigmaydp, t, p, k, ip);
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
        dJydy_Sengupta2015( dJydy,  iy,  p,  k,  y,  sigmay,  my);
    }

    virtual void fdJydy_colptrs(sunindextype *colptrs, int index) override {
        dJydy_colptrs_Sengupta2015(colptrs, index);
    }

    virtual void fdJydy_rowvals(sunindextype *rowvals, int index) override {
        dJydy_rowvals_Sengupta2015(rowvals, index);
    }


    virtual void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip) override {
        dwdp_Sengupta2015( dwdp,  t,  x,  p,  k,  h,  w,  tcl,  dtcldp,  ip);
    }


    virtual void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {
        dwdx_Sengupta2015( dwdx,  t,  x,  p,  k,  h,  w,  tcl);
    }


    virtual void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdw_Sengupta2015( dxdotdw,  t,  x,  p,  k,  h,  w);
    }


    virtual void fdxdotdw_colptrs(sunindextype *colptrs) override {
        dxdotdw_colptrs_Sengupta2015(colptrs);
    }


    virtual void fdxdotdw_rowvals(sunindextype *rowvals) override {
        dxdotdw_rowvals_Sengupta2015(rowvals);
    }


    virtual void fdxdotdp(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w) override {
        dxdotdp_Sengupta2015( dxdotdp,  t,  x,  p,  k,  h,  ip,  w);
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
        dydx_Sengupta2015(dydx, t, x, p, k, h, w, dwdx);
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
        dydp_Sengupta2015(dydp, t, x, p, k, h, ip, w, dwdp);
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
        sigmay_Sengupta2015(sigmay, t, p, k);
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
        sx0_Sengupta2015(sx0, t, x0, p, k, ip);
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
        sx0_fixedParameters_Sengupta2015(sx0, t, x0, p, k, ip);
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
        w_Sengupta2015( w,  t,  x,  p,  k,  h,  tcl);
    }


    /** model specific implementation of fx0
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0(realtype *x0, const realtype t, const realtype *p,
                     const realtype *k) override {
        x0_Sengupta2015(x0, t, p, k);
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
        x0_fixedParameters_Sengupta2015(x0, t, p, k);
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
        xdot_Sengupta2015(xdot, t, x, p, k, h, w);
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
        y_Sengupta2015(y, t, x, p, k, h, w);
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
        x_solver_Sengupta2015( x_solver,  x_rdata);
    }


    virtual void ftotal_cl(realtype *total_cl, const realtype *x_rdata) override {
        total_cl_Sengupta2015( total_cl,  x_rdata);
    }


    /**
     * @brief Get names of the model parameters
     * @return the names
     */
    virtual std::vector<std::string> getParameterNames() const override {
        return std::vector<std::string>{"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k2",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k2",
"v2",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",
"k1",
"v1",};
    }

    /**
     * @brief Get names of the model states
     * @return the names
     */
    virtual std::vector<std::string> getStateNames() const override {
        return std::vector<std::string>{"cis-Aconitate",
"car_mat",
"Xylulose5P",
"Unbranched alpha(1,4)polymer",
"UREA",
"UDP-galactose1phosphate uridylyl transferase",
"UDP-Glucose",
"UDP glucose4epimerase",
"UDP Glucose Phosphorylase",
"UDP",
"Triokinase",
"Triglyceride",
"Trehalose",
"Trehalase",
"Transketolase",
"Transaldolase",
"Thiolase",
"TPP",
"Sucrose",
"Sucrase",
"Succinate dehydrogenase",
"Succinate",
"Sedoheptulase7P",
"S CoA synthase",
"S CoA",
"Ribose5P",
"QH2",
"Q",
"Pyruvate",
"Pi",
"Pi",
"Phosphopentose isomerase",
"Phosphohexo isomerase",
"Phosphogluco mutase",
"PYRUVATE KINASE",
"PYRUVATE DEHYDROGENASE",
"PYRUVATE CARBOXYLASE",
"PPi",
"PPi",
"PHOSPHO TRIOSE ISOMERASE",
"PGA MUTASE",
"PGA KINASE",
"PEP CARBOXYKINASE",
"PEP",
"Ornithine",
"OXALOACETATE",
"OAA",
"NAPDH",
"NADPH",
"NADP+",
"NADH",
"NADH",
"NADH",
"NAD+",
"NAD+",
"Mg2+",
"Mg2+",
"Malate dehydrogenase",
"Malate",
"Limit Dextrin",
"Lactose",
"Lactonase",
"Lactate",
"LIPOATE",
"Isocitrate dehydrogenase",
"Isocitrate",
"Hexokinase",
"HMGCoA",
"HMG-CoA Synthase",
"HCO3-+ NH4+",
"H2O",
"H2O",
"H20",
"H+",
"H+",
"Glygogen Phosphorylase",
"Glycosyl-4,6-Transferase",
"Glycosyl transferase",
"Glycogen phosphorylase",
"Glycogen Synthase",
"Glycogen Primer",
"Glycogen",
"Glycerol3P",
"Glyceraldehyde",
"Glutamine",
"Glutamine",
"Glutaminase",
"Glutamate",
"Glutamate",
"Glucose",
"Galactose1P",
"Galactose",
"Galactokinase",
"GTP",
"GDP",
"GDP",
"GA3P",
"G6PDehydrogenase",
"G6P",
"G3P DEHYDROGENASE",
"G1P",
"G-6-P Phosphatase",
"Fumarate",
"Fumarase",
"Fructose",
"Fructokinase",
"FeS",
"FeS",
"FMN",
"FADH2",
"FAD",
"FA",
"F6P",
"F1P",
"F1,6P",
"F1,6BISPHOSPHATASE",
"Erythrose4P",
"Epimerase",
"Enoyl-CoA hydratase",
"ENOLASE",
"DHAP",
"D-Ribulose5P",
"D-Beta-Hydroxybutyrate",
"D Betahydroxybutyrate dehydrogenase",
"Cyb",
"Cya-a3",
"CyC2",
"CyC",
"CoA-SH",
"CoA-SH",
"Co-Ash",
"Citrullyl AMP intermediate",
"Citrulline",
"Citrate Synthase",
"Citrate",
"Carnitine_cyt",
"Carbamoyl phosphate synthetase I",
"Carbamoyl phosphate",
"CO2",
"CO2",
"CAC",
"C8Acyl-CoA",
"C8 L-3-hydroxyacyl-CoA",
"C8 Ketoacyl-CoA",
"C8 2-trans-enoyl-CoA",
"C6Acyl-CoA",
"C6 L-3-hydroxyacyl-CoA",
"C6 Ketoacyl-CoA",
"C6 2-trans-enoyl-CoA",
"C4Acyl-CoA",
"C4 L-3-hydroxyacyl-CoA",
"C4 Ketoacyl-CoA",
"C4 2-trans-enoyl-CoA",
"C22car_ims",
"C22car_ims",
"C22Acyl-CoA",
"C22 L-3-hydroxyacyl-CoA",
"C22 Ketoacyl-CoA",
"C22 AcylCoA_cyt",
"C22 2-trans-enoyl-CoA",
"C20car_ims",
"C20car_ims",
"C20Acyl-CoA",
"C20 L-3-hydroxyacyl-CoA",
"C20 Ketoacyl-CoA",
"C20 AcylCoA_cyt",
"C20 2-trans-enoyl-CoA",
"C18car_ims",
"C18car_ims",
"C18Acyl-CoA",
"C18 L-3-hydroxyacyl-CoA",
"C18 Ketoacyl-CoA",
"C18 AcylCoA_cyt",
"C18 2-trans-enoyl-CoA",
"C16car_ims",
"C16car_ims",
"C16Acyl-CoA",
"C16 L-3-hydroxyacyl-CoA",
"C16 Ketoacyl-CoA",
"C16 AcylCoA_cyt",
"C16 2-trans-enoyl-CoA",
"C14car_ims",
"C14car_ims",
"C14Acyl-CoA",
"C14 L-3-hydroxyacyl-CoA",
"C14 Ketoacyl-CoA",
"C14 AcylCoA_cyt",
"C14 2-trans-enoyl-CoA",
"C12Acyl-CoA",
"C12 L-3-hydroxyacyl-CoA",
"C12 Ketoacyl-CoA",
"C12 2-trans-enoyl-CoA",
"C10Acyl-CoA",
"C10 L-3-hydroxyacyl-CoA",
"C10 Ketoacyl-CoA",
"C10 2-trans-enoyl-CoA",
"Beta-hydroxyacyl-CoA dehydrogenase",
"Beta-KetoacylCoA dehydrogenase",
"Aspartate aminotransferase",
"Aspartate",
"Arginosuccinate",
"Arginine",
"Amino acids",
"Alpha1,6-Glycosidase",
"Alpha-KG dehydrgenase complex",
"Alpha-KG",
"Alpha keto acid",
"Alpha KG",
"Aldolase",
"Alanine",
"Acyl-CoA dehydrogenase",
"Aconitase",
"AcetylCoA",
"Acetone",
"AcetoacetylCoA",
"Acetoacetate Decarboxylase",
"Acetoacetate",
"ATP",
"ATP",
"AMP",
"ADP",
"ADP",
"A CoA",
"6PGluconate dehydrogenase",
"6PGluconate",
"6PGDL",
"3-PGA",
"2e",
"2NADH",
"2NAD+",
"2H+",
"2H+",
"2ATP",
"2ADP",
"2-PGA",
"1/2O2",
"1,3-BiPGA",
"Creatine phosphate",
"Creatine",
"Creatine Kinase",};
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
"x194",
"x195",
"x196",
"x197",
"x198",
"x199",
"x200",
"x201",
"x202",
"x203",
"x204",
"x205",
"x206",
"x207",
"x208",
"x209",
"x210",
"x211",
"x212",
"x213",
"x214",
"x215",
"x216",
"x217",
"x218",
"x219",
"x220",
"x221",
"x222",
"x223",
"x224",
"x225",
"x226",
"x227",
"x228",
"x229",
"x230",
"x231",
"x232",
"x233",
"x234",
"x235",
"x236",
"x237",
"x238",
"x239",};
    }

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getParameterIds() const override {
        return std::vector<std::string>{"re99_k1",
"re7_k1",
"re8_k1",
"re14_k1",
"re37_k1",
"re11_k1",
"re12_k1",
"re81_k1",
"re85_k1",
"re27_k1",
"re36_k1",
"re40_k1",
"re38_k1",
"re49_k1",
"re68_k1",
"re70_k1",
"re71_k1",
"re76_k1",
"re63_k1",
"re75_k1",
"re67_k1",
"re65_k1",
"re61_k1",
"re69_k1",
"re59_k1",
"re77_k1",
"re82_k1",
"re86_k1",
"re90_k1",
"re93_k1",
"re1_k2",
"re1_k1",
"re10_k1",
"re9_k1",
"re24_k1",
"re13_k1",
"re150_k1",
"re151_k1",
"re152_k1",
"re153_k1",
"re154_k1",
"re155_k1",
"re156_k1",
"re157_k1",
"re158_k1",
"re100_k1",
"re78_k1",
"re104_k1",
"re105_k1",
"re106_k1",
"re108_k1",
"re107_k1",
"re62_k1",
"re58_k1",
"re66_k1",
"re5_k1",
"re2_k1",
"re101_k1",
"re102_k1",
"re103_k1",
"re109_k1",
"re60_k1",
"re64_k1",
"re72_k1",
"re73_k1",
"re74_k1",
"re89_k1",
"re89_v1",
"re25_k1",
"re25_v1",
"re4_k1",
"re4_v1",
"re3_k1",
"re20_k1",
"re20_v1",
"re21_k1",
"re21_v1",
"re22_k1",
"re22_v1",
"re33_k1",
"re33_v1",
"re39_k1",
"re39_v1",
"re6_k1",
"re6_v1",
"re29_k1",
"re29_v1",
"re31_k1",
"re31_v1",
"re15_k1",
"re15_v1",
"re23_k1",
"re23_v1",
"re18_k1",
"re18_v1",
"re19_k1",
"re19_v1",
"re32_k1",
"re32_v1",
"re26_k1",
"re26_v1",
"re41_k1",
"re41_v1",
"re28_k1",
"re28_v1",
"re30_k1",
"re30_v1",
"re35_k1",
"re35_v1",
"re34_k1",
"re34_v1",
"re42_k1",
"re42_v1",
"re164_k1",
"re168_k1",
"re167_k1",
"re166_k1",
"re165_k1",
"re159_k1",
"re160_k1",
"re161_k1",
"re162_k1",
"re163_k1",
"re169_k1",
"re170_k1",
"re171_k1",
"re172_k1",
"re173_k1",
"re43_k1",
"re43_v1",
"re44_k1",
"re44_v1",
"re45_k1",
"re45_v1",
"re46_k1",
"re46_v1",
"re47_k1",
"re47_v1",
"re48_k1",
"re48_v1",
"re53_k1",
"re53_v1",
"re54_k1",
"re54_v1",
"re55_k1",
"re55_v1",
"re56_k1",
"re56_v1",
"re57_k1",
"re57_v1",
"re80_k1",
"re80_v1",
"re79_k1",
"re79_v1",
"re83_k1",
"re83_v1",
"re84_k1",
"re84_v1",
"re87_k1",
"re87_v1",
"re88_k1",
"re88_v1",
"re91_k1",
"re91_v1",
"re92_k1",
"re92_v1",
"re94_k1",
"re94_v1",
"re95_k1",
"re95_v1",
"re98_k1",
"re98_v1",
"re50_k1",
"re50_v1",
"re52_k1",
"re52_v1",
"re51_k1",
"re51_v1",
"re16_k1",
"re16_v1",
"re17_k1",
"re17_v1",
"re113_k1",
"re113_v1",
"re116_k1",
"re116_v1",
"re120_k1",
"re120_v1",
"re123_k1",
"re123_v1",
"re128_k1",
"re128_v1",
"re131_k1",
"re131_v1",
"re135_k1",
"re135_v1",
"re138_k2",
"re138_v2",
"re138_k1",
"re138_v1",
"re144_k1",
"re144_v1",
"re97_k1",
"re97_v1",
"re96_k1",
"re96_v1",
"re147_k1",
"re147_v1",
"re110_k1",
"re110_v1",
"re117_k1",
"re117_v1",
"re125_k1",
"re125_v1",
"re124_k1",
"re124_v1",
"re141_k1",
"re141_v1",
"re132_k1",
"re132_v1",
"re140_k1",
"re140_v1",
"re139_k1",
"re139_v1",
"re149_k1",
"re149_v1",
"re148_k1",
"re148_v1",
"re111_k1",
"re111_v1",
"re114_k1",
"re114_v1",
"re118_k1",
"re118_v1",
"re121_k1",
"re121_v1",
"re126_k1",
"re126_v1",
"re129_k1",
"re129_v1",
"re133_k1",
"re133_v1",
"re136_k1",
"re136_v1",
"re142_k1",
"re142_v1",
"re145_k1",
"re145_v1",
"re112_k1",
"re112_v1",
"re115_k1",
"re115_v1",
"re119_k1",
"re119_v1",
"re122_k1",
"re122_v1",
"re127_k1",
"re127_v1",
"re130_k1",
"re130_v1",
"re134_k1",
"re134_v1",
"re137_k1",
"re137_v1",
"re143_k1",
"re143_v1",
"re146_k1",
"re146_v1",
"re175_k1",
"re175_v1",
"re182_k1",
"re182_v1",};
    }

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    virtual std::vector<std::string> getStateIds() const override {
        return std::vector<std::string>{"s4",
"s136",
"s188",
"s253",
"s49",
"s308",
"s335",
"s13",
"s201",
"s199",
"s294",
"s292",
"s297",
"s298",
"s337",
"s248",
"s362",
"s94",
"s284",
"s283",
"s26",
"s7",
"s244",
"s58",
"s6",
"s187",
"s388",
"s267",
"s11",
"s44",
"s344",
"s193",
"s79",
"s249",
"s45",
"s91",
"s90",
"s200",
"s363",
"s36",
"s43",
"s86",
"s88",
"s10",
"s34",
"s81",
"s50",
"s190",
"s185",
"s334",
"s67",
"s351",
"s381",
"s93",
"s355",
"s333",
"s357",
"s28",
"s9",
"s252",
"s307",
"s189",
"s296",
"s95",
"s56",
"s52",
"s16",
"s322",
"s327",
"s25",
"s64",
"s347",
"s32",
"s329",
"s361",
"s255",
"s250",
"s203",
"s251",
"s336",
"s197",
"s198",
"s51",
"s293",
"s340",
"s378",
"s18",
"s341",
"s379",
"s71",
"s302",
"s306",
"s305",
"s48",
"s241",
"s352",
"s35",
"s186",
"s234",
"s85",
"s195",
"s77",
"s40",
"s27",
"s286",
"s287",
"s374",
"s377",
"s269",
"s353",
"s356",
"s53",
"s3",
"s285",
"s73",
"s31",
"s247",
"s194",
"s365",
"s87",
"s82",
"s183",
"s325",
"s330",
"s265",
"s258",
"s263",
"s262",
"s238",
"s350",
"s92",
"s38",
"s342",
"s54",
"s2",
"s126",
"s41",
"s33",
"s237",
"s349",
"s122",
"s109",
"s112",
"s111",
"s113",
"s110",
"s115",
"s116",
"s114",
"s117",
"s120",
"s119",
"s121",
"s127",
"s367",
"s15",
"s22",
"s37",
"s125",
"s19",
"s129",
"s371",
"s42",
"s68",
"s66",
"s128",
"s69",
"s131",
"s370",
"s65",
"s72",
"s80",
"s130",
"s70",
"s133",
"s369",
"s83",
"s96",
"s89",
"s132",
"s97",
"s135",
"s368",
"s84",
"s99",
"s100",
"s134",
"s98",
"s101",
"s104",
"s103",
"s105",
"s102",
"s107",
"s108",
"s106",
"s366",
"s332",
"s24",
"s23",
"s39",
"s343",
"s12",
"s254",
"s57",
"s5",
"s17",
"s21",
"s291",
"s14",
"s364",
"s354",
"s326",
"s324",
"s321",
"s328",
"s323",
"s63",
"s345",
"s47",
"s46",
"s346",
"s348",
"s192",
"s182",
"s181",
"s8",
"s259",
"s300",
"s301",
"s358",
"s389",
"s29",
"s30",
"s75",
"s256",
"s74",
"s124",
"s123",
"s400",};
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
"y194",
"y195",
"y196",
"y197",
"y198",
"y199",
"y200",
"y201",
"y202",
"y203",
"y204",
"y205",
"y206",
"y207",
"y208",
"y209",
"y210",
"y211",
"y212",
"y213",
"y214",
"y215",
"y216",
"y217",
"y218",
"y219",
"y220",
"y221",
"y222",
"y223",
"y224",
"y225",
"y226",
"y227",
"y228",
"y229",
"y230",
"y231",
"y232",
"y233",
"y234",
"y235",
"y236",
"y237",
"y238",
"y239",};
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
