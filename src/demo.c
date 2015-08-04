#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <math.h>
#include "utilities.h"

/* define useful macros */
static double* parms;
static int nAgeGroups;
static int years;
static int reverse;

#define BIRTHS(i) (parms[i])
#define DEATHS(i) (parms[years + i])
#define AGEING(i) (parms[2 * years + i])

#define Ndot(i) (ydot[i])
#define N(i) (y[i])
#define totalN (yout[0])

/* initialiser -- assign R variables to C variables */
void demo_initmod(void (* odeparms)(int *, double *))
{
    DL_FUNC get_deSolve_gparms;
    get_deSolve_gparms = R_GetCCallable("deSolve","get_deSolve_gparms");
    SEXP gparms = get_deSolve_gparms();

    SEXP Rbirths = getListElement(gparms, "births");
    if (Rbirths == NULL) {
        error("No births given");
    }
    int *births = INTEGER(Rbirths);
    SEXP Rdeaths = getListElement(gparms, "deaths");
    if (Rdeaths == NULL) {
        error("No deaths given");
    }
    int *deaths = INTEGER(Rdeaths);
    years = LENGTH(Rbirths);
    SEXP Ragewidths = getListElement(gparms, "agewidths");
    double* ageing;
    double* agewidths;
    int age_widths_not_rates = 0;
    if (Ragewidths == NULL) {
        SEXP Rageing = getListElement(gparms, "ageing");
        if (Rageing == NULL) {
            error("No age widths given");
        }
        nAgeGroups = LENGTH(Rageing) + 1;
        ageing = REAL(Rageing);
    } else {
        nAgeGroups = LENGTH(Ragewidths) + 1;
        agewidths = REAL(Ragewidths);
        age_widths_not_rates = 1;
    }

    SEXP Rreverse = getListElement(gparms, "reverse");
    if (Rreverse == NULL) {
        reverse = 0;
    } else {
        reverse = LOGICAL(Rreverse)[0];
    }
    int N = 2 * years + nAgeGroups - 1;
    parms = malloc(N * sizeof(double));


    /* printf("years:%i nAgeGroups:%i\n", years, nAgeGroups); */

    for (size_t i = 0; i < years; ++i) {
        /* printf("i:%lu births:%f\n", i, births[i]); */
        BIRTHS(i) = births[(reverse ? years - i - 1 : i)];
    }
    for (size_t i = 0; i < years; ++i) {
        /* printf("i:%lu deaths:%f\n", i, deaths[i]); */
        DEATHS(i) = deaths[(reverse ? years - i - 1 : i)];
    }
    for (size_t i = 0; i < (nAgeGroups - 1); ++i) {
        /* printf("i:%lu agewidths:%f\n", i, agewidths[i]); */
        if (age_widths_not_rates) {
            AGEING(i) = 1 / (agewidths[i] + .0);
        } else {
            AGEING(i) = ageing[i];
        }
    }
}

/* derivatives and 1 output variable */
void demo_derivs (int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip)
{
    if (ip[0] < 1) error("nout should be at least 1");

    size_t year = (*t) / 1;

    for (size_t i = 0; i < nAgeGroups; ++i) {
        Ndot(i) = (i == 0 ? BIRTHS(year) : 0) -
            (i < (nAgeGroups - 1) ? AGEING(i) * N(i) : 0) +
            (i > 0 ? AGEING(i - 1) * N(i - 1) : 0) -
            (i == (nAgeGroups - 1) ? DEATHS(year) : 0);
        if (reverse) {
            Ndot(i) = -Ndot(i);
        }
        /* printf("i:%lu t:%f year:%lu births:%f ageing:%f deaths:%f N:%f change:%f totalN:%f\n", i, *t, year, BIRTHS(year), (i < (nAgeGroups - 1) ? AGEING(i) : 0), DEATHS(year), N(i), Ndot(i), totalN); */
    }

    totalN = 0;
    for (size_t i = 0; i < nAgeGroups; ++i) {
        totalN = totalN + N(i);
    }

}
