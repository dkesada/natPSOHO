// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// create_causlist_cpp
Rcpp::List create_causlist_cpp(Rcpp::List& net, unsigned int size, StringVector& ordering);
RcppExport SEXP _psoho_create_causlist_cpp(SEXP netSEXP, SEXP sizeSEXP, SEXP orderingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type net(netSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< StringVector& >::type ordering(orderingSEXP);
    rcpp_result_gen = Rcpp::wrap(create_causlist_cpp(net, size, ordering));
    return rcpp_result_gen;
END_RCPP
}
// cl_to_arc_matrix_cpp
Rcpp::CharacterMatrix cl_to_arc_matrix_cpp(Rcpp::List& cl, Rcpp::CharacterVector& ordering, Rcpp::NumericMatrix& counters, unsigned int rows);
RcppExport SEXP _psoho_cl_to_arc_matrix_cpp(SEXP clSEXP, SEXP orderingSEXP, SEXP countersSEXP, SEXP rowsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type cl(clSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector& >::type ordering(orderingSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type counters(countersSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type rows(rowsSEXP);
    rcpp_result_gen = Rcpp::wrap(cl_to_arc_matrix_cpp(cl, ordering, counters, rows));
    return rcpp_result_gen;
END_RCPP
}
// rename_nodes_cpp
Rcpp::StringVector rename_nodes_cpp(Rcpp::StringVector& nodes, unsigned int size);
RcppExport SEXP _psoho_rename_nodes_cpp(SEXP nodesSEXP, SEXP sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector& >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type size(sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(rename_nodes_cpp(nodes, size));
    return rcpp_result_gen;
END_RCPP
}
// initialize_vl_cpp
Rcpp::List initialize_vl_cpp(StringVector& ordering, unsigned int size);
RcppExport SEXP _psoho_initialize_vl_cpp(SEXP orderingSEXP, SEXP sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector& >::type ordering(orderingSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type size(sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(initialize_vl_cpp(ordering, size));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_psoho_create_causlist_cpp", (DL_FUNC) &_psoho_create_causlist_cpp, 3},
    {"_psoho_cl_to_arc_matrix_cpp", (DL_FUNC) &_psoho_cl_to_arc_matrix_cpp, 4},
    {"_psoho_rename_nodes_cpp", (DL_FUNC) &_psoho_rename_nodes_cpp, 2},
    {"_psoho_initialize_vl_cpp", (DL_FUNC) &_psoho_initialize_vl_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_psoho(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
