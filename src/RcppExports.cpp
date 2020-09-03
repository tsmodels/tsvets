// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// vets_cpp_llh_equicor
Rcpp::List vets_cpp_llh_equicor(NumericVector model, arma::mat Amat, arma::sp_mat Fmat, arma::sp_mat Hmat, arma::sp_mat Gmat, arma::mat States, arma::mat Y, arma::mat X, arma::mat beta);
RcppExport SEXP _tsvets_vets_cpp_llh_equicor(SEXP modelSEXP, SEXP AmatSEXP, SEXP FmatSEXP, SEXP HmatSEXP, SEXP GmatSEXP, SEXP StatesSEXP, SEXP YSEXP, SEXP XSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Amat(AmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Fmat(FmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Hmat(HmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Gmat(GmatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type States(StatesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(vets_cpp_llh_equicor(model, Amat, Fmat, Hmat, Gmat, States, Y, X, beta));
    return rcpp_result_gen;
END_RCPP
}
// vets_cpp_llh_diagonal
Rcpp::List vets_cpp_llh_diagonal(NumericVector model, arma::mat Amat, arma::sp_mat Fmat, arma::sp_mat Hmat, arma::sp_mat Gmat, arma::mat States, arma::mat Y, arma::mat X, arma::mat beta);
RcppExport SEXP _tsvets_vets_cpp_llh_diagonal(SEXP modelSEXP, SEXP AmatSEXP, SEXP FmatSEXP, SEXP HmatSEXP, SEXP GmatSEXP, SEXP StatesSEXP, SEXP YSEXP, SEXP XSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Amat(AmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Fmat(FmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Hmat(HmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Gmat(GmatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type States(StatesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(vets_cpp_llh_diagonal(model, Amat, Fmat, Hmat, Gmat, States, Y, X, beta));
    return rcpp_result_gen;
END_RCPP
}
// vets_cpp_llh_full
Rcpp::List vets_cpp_llh_full(NumericVector model, arma::mat Amat, arma::sp_mat Fmat, arma::sp_mat Hmat, arma::sp_mat Gmat, arma::mat States, arma::mat Y, arma::mat X, arma::mat beta);
RcppExport SEXP _tsvets_vets_cpp_llh_full(SEXP modelSEXP, SEXP AmatSEXP, SEXP FmatSEXP, SEXP HmatSEXP, SEXP GmatSEXP, SEXP StatesSEXP, SEXP YSEXP, SEXP XSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Amat(AmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Fmat(FmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Hmat(HmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Gmat(GmatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type States(StatesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(vets_cpp_llh_full(model, Amat, Fmat, Hmat, Gmat, States, Y, X, beta));
    return rcpp_result_gen;
END_RCPP
}
// vets_cpp_llh_shrink
Rcpp::List vets_cpp_llh_shrink(NumericVector model, arma::mat Amat, arma::sp_mat Fmat, arma::sp_mat Hmat, arma::sp_mat Gmat, arma::mat States, arma::mat Y, arma::mat X, arma::mat beta);
RcppExport SEXP _tsvets_vets_cpp_llh_shrink(SEXP modelSEXP, SEXP AmatSEXP, SEXP FmatSEXP, SEXP HmatSEXP, SEXP GmatSEXP, SEXP StatesSEXP, SEXP YSEXP, SEXP XSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Amat(AmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Fmat(FmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Hmat(HmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Gmat(GmatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type States(StatesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(vets_cpp_llh_shrink(model, Amat, Fmat, Hmat, Gmat, States, Y, X, beta));
    return rcpp_result_gen;
END_RCPP
}
// vets_cpp_filter
Rcpp::List vets_cpp_filter(NumericVector model, arma::mat Amat, arma::sp_mat Fmat, arma::sp_mat Hmat, arma::sp_mat Gmat, arma::mat States, arma::mat Y, arma::mat X, arma::mat beta);
RcppExport SEXP _tsvets_vets_cpp_filter(SEXP modelSEXP, SEXP AmatSEXP, SEXP FmatSEXP, SEXP HmatSEXP, SEXP GmatSEXP, SEXP StatesSEXP, SEXP YSEXP, SEXP XSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Amat(AmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Fmat(FmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Hmat(HmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Gmat(GmatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type States(StatesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(vets_cpp_filter(model, Amat, Fmat, Hmat, Gmat, States, Y, X, beta));
    return rcpp_result_gen;
END_RCPP
}
// vets_cpp_predict
Rcpp::List vets_cpp_predict(Rcpp::NumericVector model, arma::mat S, arma::mat Amat, arma::sp_mat Fmat, arma::sp_mat Hmat, arma::sp_mat Gmat, arma::vec istate, arma::mat X, arma::mat beta);
RcppExport SEXP _tsvets_vets_cpp_predict(SEXP modelSEXP, SEXP SSEXP, SEXP AmatSEXP, SEXP FmatSEXP, SEXP HmatSEXP, SEXP GmatSEXP, SEXP istateSEXP, SEXP XSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Amat(AmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Fmat(FmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Hmat(HmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Gmat(GmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type istate(istateSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(vets_cpp_predict(model, S, Amat, Fmat, Hmat, Gmat, istate, X, beta));
    return rcpp_result_gen;
END_RCPP
}
// vets_cpp_simulate
Rcpp::List vets_cpp_simulate(Rcpp::NumericVector model, arma::cube E, arma::mat Amat, arma::sp_mat Fmat, arma::sp_mat Hmat, arma::sp_mat Gmat, arma::vec istate, arma::mat X, arma::mat beta);
RcppExport SEXP _tsvets_vets_cpp_simulate(SEXP modelSEXP, SEXP ESEXP, SEXP AmatSEXP, SEXP FmatSEXP, SEXP HmatSEXP, SEXP GmatSEXP, SEXP istateSEXP, SEXP XSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type E(ESEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Amat(AmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Fmat(FmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Hmat(HmatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Gmat(GmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type istate(istateSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(vets_cpp_simulate(model, E, Amat, Fmat, Hmat, Gmat, istate, X, beta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tsvets_vets_cpp_llh_equicor", (DL_FUNC) &_tsvets_vets_cpp_llh_equicor, 9},
    {"_tsvets_vets_cpp_llh_diagonal", (DL_FUNC) &_tsvets_vets_cpp_llh_diagonal, 9},
    {"_tsvets_vets_cpp_llh_full", (DL_FUNC) &_tsvets_vets_cpp_llh_full, 9},
    {"_tsvets_vets_cpp_llh_shrink", (DL_FUNC) &_tsvets_vets_cpp_llh_shrink, 9},
    {"_tsvets_vets_cpp_filter", (DL_FUNC) &_tsvets_vets_cpp_filter, 9},
    {"_tsvets_vets_cpp_predict", (DL_FUNC) &_tsvets_vets_cpp_predict, 9},
    {"_tsvets_vets_cpp_simulate", (DL_FUNC) &_tsvets_vets_cpp_simulate, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_tsvets(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
