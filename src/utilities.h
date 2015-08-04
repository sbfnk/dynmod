#ifndef UTILITIES_H
#define UTILITIES_H

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <math.h>

SEXP getListElement(SEXP list, const char *str);

#endif
