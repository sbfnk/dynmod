#include "utilities.h"

/* get an element of a list, return null if it doesn't exist */
SEXP getListElement(SEXP list, const char *str) {
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    int i;
    for (i = 0; i < length(list); i++)
        if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
            elmt = VECTOR_ELT(list, i);
            break;
        }

    if (elmt == R_NilValue)
        return NULL;
    else
        return elmt;
}
