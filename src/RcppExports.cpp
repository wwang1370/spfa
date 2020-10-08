// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// spfa_main
Rcpp::List spfa_main(const arma::mat& shortpar, arma::mat& dat, const arma::uvec& type, arma::uword n_basis, double lmbd, arma::uword n_quad, arma::uword maxit_em, arma::uword maxit_mstep, arma::uword maxit_start, double tol_em, double tol_mstep, double tol_start);
RcppExport SEXP _spfa_spfa_main(SEXP shortparSEXP, SEXP datSEXP, SEXP typeSEXP, SEXP n_basisSEXP, SEXP lmbdSEXP, SEXP n_quadSEXP, SEXP maxit_emSEXP, SEXP maxit_mstepSEXP, SEXP maxit_startSEXP, SEXP tol_emSEXP, SEXP tol_mstepSEXP, SEXP tol_startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type shortpar(shortparSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type dat(datSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type type(typeSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_basis(n_basisSEXP);
    Rcpp::traits::input_parameter< double >::type lmbd(lmbdSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_quad(n_quadSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type maxit_em(maxit_emSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type maxit_mstep(maxit_mstepSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type maxit_start(maxit_startSEXP);
    Rcpp::traits::input_parameter< double >::type tol_em(tol_emSEXP);
    Rcpp::traits::input_parameter< double >::type tol_mstep(tol_mstepSEXP);
    Rcpp::traits::input_parameter< double >::type tol_start(tol_startSEXP);
    rcpp_result_gen = Rcpp::wrap(spfa_main(shortpar, dat, type, n_basis, lmbd, n_quad, maxit_em, maxit_mstep, maxit_start, tol_em, tol_mstep, tol_start));
    return rcpp_result_gen;
END_RCPP
}
// spfa_score
arma::vec spfa_score(const arma::mat& shortpar, arma::mat& dat, arma::uword n_basis, arma::uword n_quad);
RcppExport SEXP _spfa_spfa_score(SEXP shortparSEXP, SEXP datSEXP, SEXP n_basisSEXP, SEXP n_quadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type shortpar(shortparSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type dat(datSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_basis(n_basisSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_quad(n_quadSEXP);
    rcpp_result_gen = Rcpp::wrap(spfa_score(shortpar, dat, n_basis, n_quad));
    return rcpp_result_gen;
END_RCPP
}
// marg_lik
arma::vec marg_lik(const arma::mat& shortpar, arma::mat& dat, arma::uword n_basis, arma::uword n_quad);
RcppExport SEXP _spfa_marg_lik(SEXP shortparSEXP, SEXP datSEXP, SEXP n_basisSEXP, SEXP n_quadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type shortpar(shortparSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type dat(datSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_basis(n_basisSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_quad(n_quadSEXP);
    rcpp_result_gen = Rcpp::wrap(marg_lik(shortpar, dat, n_basis, n_quad));
    return rcpp_result_gen;
END_RCPP
}
// xv_risk
arma::vec xv_risk(const arma::mat& shortpar, arma::mat& dat, arma::uword n_basis, arma::uword n_quad, arma::uword order);
RcppExport SEXP _spfa_xv_risk(SEXP shortparSEXP, SEXP datSEXP, SEXP n_basisSEXP, SEXP n_quadSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type shortpar(shortparSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type dat(datSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_basis(n_basisSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_quad(n_quadSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(xv_risk(shortpar, dat, n_basis, n_quad, order));
    return rcpp_result_gen;
END_RCPP
}
// cond_dns
arma::mat cond_dns(arma::vec shortpar, arma::vec y, arma::vec x, arma::uword n_basis, arma::uword n_quad);
RcppExport SEXP _spfa_cond_dns(SEXP shortparSEXP, SEXP ySEXP, SEXP xSEXP, SEXP n_basisSEXP, SEXP n_quadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type shortpar(shortparSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_basis(n_basisSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_quad(n_quadSEXP);
    rcpp_result_gen = Rcpp::wrap(cond_dns(shortpar, y, x, n_basis, n_quad));
    return rcpp_result_gen;
END_RCPP
}
// cubic_bspl
arma::mat cubic_bspl(arma::vec x, arma::uword n_basis);
RcppExport SEXP _spfa_cubic_bspl(SEXP xSEXP, SEXP n_basisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_basis(n_basisSEXP);
    rcpp_result_gen = Rcpp::wrap(cubic_bspl(x, n_basis));
    return rcpp_result_gen;
END_RCPP
}
// cubic_bspl0
arma::mat cubic_bspl0(arma::vec x, arma::uword n_basis, double x0);
RcppExport SEXP _spfa_cubic_bspl0(SEXP xSEXP, SEXP n_basisSEXP, SEXP x0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_basis(n_basisSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    rcpp_result_gen = Rcpp::wrap(cubic_bspl0(x, n_basis, x0));
    return rcpp_result_gen;
END_RCPP
}
// gl_quad
Rcpp::List gl_quad(arma::uword n_quad, double lower, double upper);
RcppExport SEXP _spfa_gl_quad(SEXP n_quadSEXP, SEXP lowerSEXP, SEXP upperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type n_quad(n_quadSEXP);
    Rcpp::traits::input_parameter< double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< double >::type upper(upperSEXP);
    rcpp_result_gen = Rcpp::wrap(gl_quad(n_quad, lower, upper));
    return rcpp_result_gen;
END_RCPP
}
// reduce_par
arma::vec reduce_par(arma::uword n_basis, arma::vec par);
RcppExport SEXP _spfa_reduce_par(SEXP n_basisSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type n_basis(n_basisSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(reduce_par(n_basis, par));
    return rcpp_result_gen;
END_RCPP
}
// extend_par
arma::vec extend_par(arma::uword n_basis, arma::vec shortpar);
RcppExport SEXP _spfa_extend_par(SEXP n_basisSEXP, SEXP shortparSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type n_basis(n_basisSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type shortpar(shortparSEXP);
    rcpp_result_gen = Rcpp::wrap(extend_par(n_basis, shortpar));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spfa_spfa_main", (DL_FUNC) &_spfa_spfa_main, 12},
    {"_spfa_spfa_score", (DL_FUNC) &_spfa_spfa_score, 4},
    {"_spfa_marg_lik", (DL_FUNC) &_spfa_marg_lik, 4},
    {"_spfa_xv_risk", (DL_FUNC) &_spfa_xv_risk, 5},
    {"_spfa_cond_dns", (DL_FUNC) &_spfa_cond_dns, 5},
    {"_spfa_cubic_bspl", (DL_FUNC) &_spfa_cubic_bspl, 2},
    {"_spfa_cubic_bspl0", (DL_FUNC) &_spfa_cubic_bspl0, 3},
    {"_spfa_gl_quad", (DL_FUNC) &_spfa_gl_quad, 3},
    {"_spfa_reduce_par", (DL_FUNC) &_spfa_reduce_par, 2},
    {"_spfa_extend_par", (DL_FUNC) &_spfa_extend_par, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_spfa(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
