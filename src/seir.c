// todo: do it all by transition
// inconsitencies: boosted_immunity, boosterYears etc; antibody dvelopment & uptake

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <math.h>
#include "utilities.h"

/* define useful macros */
static int nAgeGroups;
static int nE;
static int nC;
static int nI;
static int termtime;
static int sinusoidal;
static int birthTimes;
static int deathTimes;
static int bedTimes;
static int childMixingTimes;
static int uptakeTimes;
static int campaignTimes;
static int boosterCampaignTimes;
static int verbose;
static double beta;
static double pgamma;
#define gamma pgamma
static double delta;
static double epsilon;
static double rho;
static double r;
static double asymp;
static double termtime_holiday_start;
static double termtime_holiday_end;
static double sinusoidal_amplitude;
static double sinusoidal_t0;
static double sinusoidal_period;
static double birthrate;
static double deathrate;
static double vaccine_efficacy;
static double maternal_immunity;
static double antibody_development;
static double leaky_immunity;
static double boosted_immunity;
static double vaccine_waning;
static double boosted_waning;
static double hospital_entry;
static int interpolate;

static double *ageing;
static double *pd;

static double *mixing_matrix;
static double *termtime_holiday_mixing_matrix;
static double *uptake;
static double *campaign;
static double *boosterCampaign;

static int *births;
static int *deaths;
static int *beds;
static double* childMixing;

static int childMixingGroup;

static int start_step;
static int started;

int get_ivalue(int* v, double time, unsigned int steps)
{
    unsigned int current_step;

    if (interpolate)
    {
        current_step = floor(time + 0.5) - start_step;
    } else
    {
        current_step = floor(time) - start_step;
    }
    double res;

    if (steps > 0)
    {
        if (interpolate)
        {
            if (current_step < (steps - 1))
            {
                res = v[current_step + 1] *
                    (time - start_step - current_step + 0.5) -
                    v[current_step] *
                    (time - current_step - start_step - 0.5);
            } else
            {
                res = v[steps - 1];
            }
        } else
        {
            res = v[(current_step < steps) ? (current_step) : (steps - 1)];
        }
    } else
    {
        res = 0;
    }
    return(round(res));
}

double get_fvalue(double* v, double time, unsigned int steps)
{
    unsigned int current_step;

    if (interpolate)
    {
        current_step = floor(time + 0.5) - start_step;
    } else
    {
        current_step = floor(time) - start_step;
    }

    double res;

    if (steps > 0)
    {
        if (interpolate)
        {
            if (current_step < (steps - 1))
            {
                res = v[current_step + 1] *
                    (time - start_step - current_step + 0.5) -
                    v[current_step] *
                    (time - current_step - start_step - 0.5);
            } else
            {
                res = v[steps - 1];
            }
        } else
        {
            res = v[(current_step < steps) ? (current_step) : (steps - 1)];
        }
    } else
    {
        res = 0;
    }

    return(res);
}

double get_fvalues(double* v, double time, unsigned int group,
                   unsigned int steps, unsigned int groups)
{
    unsigned int current_step;

    if (interpolate)
    {
        // if interpolating, the current step starts 0.5 before the time
        current_step = floor(time + 0.5) - start_step;
    } else
    {
        // if not interpolating, the current step starts at 0
        current_step = floor(time) - start_step;
    }

    double res;

    if (steps > 0)
    {
        if (interpolate)
        {
            if (current_step < (start_step + steps - 1))
            {
                res = v[(current_step + 1) + steps * group] *
                    (time - start_step - current_step + 0.5) -
                    v[current_step + steps * group] *
                    (time - start_step - current_step - 0.5);
            } else
            {
                res = v[(steps - 1) + steps * group];
            }
        } else
        {
            // if not interpolating, we check if we have enough
            // information -- if yes, get the right matrix element
            res = v[(current_step < steps) ? (current_step + steps * group) :
                    ((steps - 1) + steps * group)];
        }
    } else
    {
        res = 0;
    }

    return(res);
}

#define NOR_MIX(i, j) (mixing_matrix[i + nAgeGroups * j])
#define HOL_MIX(i, j) (termtime_holiday_mixing_matrix[i + nAgeGroups * j])

#define B (maternal_immunity > 0 ? y[0] : 0)
#define S(i) (y[(maternal_immunity > 0 ? 1 : 0) + i])
#define E(i, j) (y[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * i + j])
#define C(i, j) (y[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * nAgeGroups + nC * i + j])
#define I(i, j) (y[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * nAgeGroups + nC * nAgeGroups + nI * i + j])
#define R(i) (y[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * nAgeGroups + nC * nAgeGroups + nI * nAgeGroups + i])
#define A(i) (y[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * nAgeGroups + nC * nAgeGroups + nI * nAgeGroups + nAgeGroups + i])
#define V(i) (y[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * nAgeGroups + nC * nAgeGroups + nI * nAgeGroups + nAgeGroups + (antibody_development > 0 ? 1 : 0) * nAgeGroups + i])
#define W(i) (y[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * nAgeGroups + nC * nAgeGroups + nI * nAgeGroups + nAgeGroups + (antibody_development > 0 ? 1 : 0) * nAgeGroups + (uptakeTimes + campaignTimes > 0 ? 1 : 0) * nAgeGroups + i])
#define H(i) (y[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * nAgeGroups + nC * nAgeGroups + nI * nAgeGroups + nAgeGroups + (antibody_development > 0 ? 1 : 0) * nAgeGroups + (uptakeTimes + campaignTimes > 0 ? 1 : 0) * nAgeGroups + (boosterCampaignTimes > 0 ? 1 : 0) * nAgeGroups + i])
#define Z(i) (y[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * nAgeGroups + nC * nAgeGroups + nI * nAgeGroups + nAgeGroups + (antibody_development > 0 ? 1 : 0) * nAgeGroups + ((uptakeTimes + campaignTimes > 0) ? 1 : 0) * nAgeGroups + (boosterCampaignTimes > 0 ? 1 : 0) * nAgeGroups + (bedTimes > 0 ? 1 : 0) * nAgeGroups + i])

#define Bdot (ydot[0])
#define Sdot(i) (ydot[(maternal_immunity > 0 ? 1 : 0) + i])
#define Edot(i, j) (ydot[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * i + j])
#define Cdot(i, j) (ydot[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * nAgeGroups + nC * i + j])
#define Idot(i, j) (ydot[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * nAgeGroups + nC * nAgeGroups + nI * i + j])
#define Rdot(i) (ydot[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * nAgeGroups + nC * nAgeGroups + nI * nAgeGroups + i])
#define Adot(i) (ydot[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * nAgeGroups + nC * nAgeGroups + nI * nAgeGroups + nAgeGroups + i])
#define Vdot(i) (ydot[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * nAgeGroups + nC * nAgeGroups + nI * nAgeGroups + nAgeGroups + (antibody_development > 0 ? 1 : 0) * nAgeGroups + i])
#define Wdot(i) (ydot[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * nAgeGroups + nC * nAgeGroups + nI * nAgeGroups + nAgeGroups + (antibody_development > 0 ? 1 : 0) * nAgeGroups + (uptakeTimes + campaignTimes > 0 ? 1 : 0) * nAgeGroups + i])
#define Zdot(i) (ydot[(maternal_immunity > 0 ? 1 : 0) + nAgeGroups + nE * nAgeGroups + nC * nAgeGroups + nI * nAgeGroups + nAgeGroups + (antibody_development > 0 ? 1 : 0) * nAgeGroups + ((uptakeTimes + campaignTimes > 0) ? 1 : 0) * nAgeGroups + (boosterCampaignTimes > 0 ? 1 : 0) * nAgeGroups + i])

#define N(i) (yout[i])
#define totalN (yout[nAgeGroups])

/* initialiser -- assign R variables to C variables */
void seir_initmod(void (* odeparms)(int *, double *))
{
    DL_FUNC get_deSolve_gparms;
    get_deSolve_gparms = R_GetCCallable("deSolve","get_deSolve_gparms");
    SEXP gparms = get_deSolve_gparms();

    SEXP Rverbose = getListElement(gparms, "verbose");
    if (!(Rverbose == NULL))
    {
        verbose = INTEGER(Rverbose)[0];
    }

    SEXP Rmixing_matrix = getListElement(gparms, "mixing");
    SEXP Rageing = getListElement(gparms, "ageing");
    if (Rmixing_matrix == NULL)
    {
        if (Rageing == NULL)
        {
            nAgeGroups = 1;
        } else
        {
            nAgeGroups = LENGTH(Rageing) + 1;
        }
    } else
    {
        SEXP mixing_dim = getAttrib( Rmixing_matrix, R_DimSymbol );
        int mixing_nrow = INTEGER(mixing_dim)[0];
        nAgeGroups = mixing_nrow;
    }

    if (Rageing == NULL)
    {
        ageing = malloc((nAgeGroups - 1) * sizeof(double));
        for (unsigned int i = 0; i < (nAgeGroups - 1); ++i)
        {
            ageing[i] = 0;
        }
    } else
    {
        ageing = REAL(Rageing);
    }


    if (Rmixing_matrix == NULL)
    {
        mixing_matrix = malloc(sizeof(double) * nAgeGroups * nAgeGroups);
        for (unsigned int i = 0; i < nAgeGroups * nAgeGroups; ++i)
        {
            mixing_matrix[i] = 1;
        }
    } else
    {
        mixing_matrix = REAL(Rmixing_matrix);
    }

    SEXP RtermtimeForcing = getListElement(gparms, "termtime.forcing");

    if (RtermtimeForcing == NULL)
    {
        termtime = 0;
    } else
    {
        termtime = 1;
    }

    SEXP RsinusoidalForcing = getListElement(gparms, "sinusoidal.forcing");

    if (RsinusoidalForcing == NULL)
    {
        sinusoidal = 0;
    } else
    {
        sinusoidal = 1;
    }

    SEXP Rbirths = getListElement(gparms, "births");
    SEXP Rmu = getListElement(gparms, "mu");
    if (Rbirths == NULL)
    {
        if (Rmu == NULL)
        {
            birthrate = 0;
        } else
        {
            birthrate = REAL(Rmu)[0];
            if (verbose)
            {
                printf("birthrate:%f\n", birthrate);
            }
        }
        birthTimes = 0;
    } else
    {
        births = INTEGER(Rbirths);
        birthTimes = LENGTH(Rbirths);
    }

    SEXP Rdeaths = getListElement(gparms, "deaths");
    if (Rdeaths == NULL)
    {
        SEXP Rmu = getListElement(gparms, "mu");
        if (Rmu == NULL)
        {
            deathrate = 0;
        } else
        {
            deathrate = REAL(Rmu)[0];
            if (verbose)
            {
                printf("deathrate:%f\n", deathrate);
            }
        }
        deathTimes = 0;
    } else
    {
        deaths = INTEGER(Rdeaths);
        deathTimes = LENGTH(Rdeaths);
    }

    SEXP Rpd = getListElement(gparms, "pd");
    if (Rpd == NULL)
    {
        pd = malloc(nAgeGroups * sizeof(double));
        for (unsigned int i = 0; i < nAgeGroups; ++i)
        {
            pd[i] = 1;
        }
    } else
    {
        pd = REAL(Rpd);
    }

    SEXP RchildMixing = getListElement(gparms, "child.mixing");
    if (RchildMixing == NULL)
    {
        childMixing = malloc(sizeof(double));
        childMixing[0] = 1;
        childMixingTimes = 1;
    } else if (LENGTH(RchildMixing) > 0)
    {
        childMixing = REAL(RchildMixing);
        childMixingTimes = LENGTH(RchildMixing);
    }

    SEXP RchildMixingGroup = getListElement(gparms, "child.mixing.group");
    if (RchildMixingGroup == NULL)
    {
        childMixingGroup = 0;
    } else
    {
        childMixingGroup = INTEGER(RchildMixingGroup)[0] - 1;
    }

    SEXP Rbeta = getListElement(gparms, "beta");
    if (Rbeta == NULL)
    {
        beta = 0;
    } else
    {
        beta = REAL(Rbeta)[0];
    }
    SEXP Rgamma = getListElement(gparms, "gamma");
    if (Rgamma == NULL)
    {
        gamma = 0;
    } else
    {
        gamma = REAL(Rgamma)[0];
    }
    SEXP Rdelta = getListElement(gparms, "delta");
    if (Rdelta == NULL)
    {
        delta = 0;
    } else
    {
        delta = REAL(Rdelta)[0];
    }

    SEXP RnE = getListElement(gparms, "nE");
    if (RnE == NULL)
    {
        error("nE must be provided");
    } else
    {
        nE = INTEGER(RnE)[0];

    }

    SEXP RnC = getListElement(gparms, "nC");
    if (RnC == NULL)
    {
        error("nC must be provided");
    } else
    {
        nC = INTEGER(RnC)[0];
    }

    SEXP RnI = getListElement(gparms, "nI");
    if (RnI == NULL)
    {
        error("nI must be provided");
    } else
    {
        nI = INTEGER(RnI)[0];
    }

    SEXP Repsilon = getListElement(gparms, "epsilon");
    if (Repsilon == NULL)
    {
        epsilon = 0;
    } else
    {
        if (nC == 0)
        {
            error("epsilon provided but nC = 0");
        } else
        {
            epsilon = REAL(Repsilon)[0];
        }
    }

    SEXP Rrho = getListElement(gparms, "rho");
    if (Rrho == NULL)
    {
        rho = 0;
    } else
    {
        if (nE == 0)
        {
            error("rho provided but nE = 0");
        } else
        {
            rho = REAL(Rrho)[0];
        }
    }

    SEXP Rasymp = getListElement(gparms, "asymp");
    if (Rasymp == NULL)
    {
        asymp = 0;
    } else
    {
        asymp = REAL(Rasymp)[0];
        if (verbose)
        {
            printf("asymp:%f\n", asymp);
        }
    }

    SEXP Rreservoir = getListElement(gparms, "reservoir");
    if (Rreservoir == NULL)
    {
        r = 0;
    } else
    {
        r = REAL(getListElement(gparms, "reservoir"))[0];
    }

    SEXP Refficacy = getListElement(gparms, "vaccine.efficacy");
    if (Refficacy == NULL)
    {
        vaccine_efficacy = 1;
    } else
    {
        vaccine_efficacy = REAL(Refficacy)[0];
    }

    SEXP Rmaternal_immunity = getListElement(gparms, "maternal.immunity");
    if (Rmaternal_immunity == NULL)
    {
        maternal_immunity = 0;
    } else
    {
        maternal_immunity = REAL(Rmaternal_immunity)[0];
    }

    SEXP Rantibody_development = getListElement(gparms, "antibody.development");
    if (Rantibody_development == NULL)
    {
        antibody_development = 0;
    } else
    {
        antibody_development = REAL(Rantibody_development)[0];
    }

    SEXP Rleaky_immunity = getListElement(gparms, "leaky.immunity");
    if (Rleaky_immunity == NULL)
    {
        leaky_immunity = 1;
    } else
    {
        leaky_immunity = REAL(Rleaky_immunity)[0];
    }

    SEXP Rboosted_immunity = getListElement(gparms, "boosted.immunity");
    if (Rboosted_immunity == NULL)
    {
        boosted_immunity = 1;
    } else
    {
        boosted_immunity = REAL(Rboosted_immunity)[0];
    }

    SEXP Rvaccine_waning = getListElement(gparms, "vaccine.waning");
    if (Rvaccine_waning == NULL)
    {
        vaccine_waning = 0;
    } else
    {
        vaccine_waning = REAL(Rvaccine_waning)[0];
    }

    SEXP Rboosted_waning = getListElement(gparms, "boosted.waning");
    if (Rboosted_waning == NULL)
    {
        boosted_waning = 0;
    } else
    {
        boosted_waning = REAL(Rboosted_waning)[0];
    }

    SEXP Rhospital_entry = getListElement(gparms, "hospital.entry");
    if (Rhospital_entry == NULL)
    {
        hospital_entry = 0;
    } else
    {
        hospital_entry = REAL(Rhospital_entry)[0];
    }

    SEXP Rinterpolate = getListElement(gparms, "interpolate");
    if (!(Rinterpolate == NULL))
    {
        interpolate = LOGICAL(Rinterpolate)[0];
    }

    if (termtime)
    {
        termtime_holiday_start = REAL(getListElement(RtermtimeForcing,
                                                     "holiday.start"))[0];
        termtime_holiday_end = REAL(getListElement(RtermtimeForcing,
                                                   "holiday.end"))[0];
        SEXP Rholiday_mixing_matrix = getListElement(RtermtimeForcing,
                                                     "holiday.mixing");
        termtime_holiday_mixing_matrix = REAL(Rholiday_mixing_matrix);
    } else
    {
        termtime = 0;
    }

    if (getListElement(gparms, "sinusoidal.forcing") != NULL)
    {
        sinusoidal = 1;
        SEXP RsinusoidalForcing = getListElement(gparms, "sinusoidal.forcing");
        sinusoidal_amplitude = REAL(getListElement(RsinusoidalForcing,
                                                   "amplitude"))[0];
        sinusoidal_t0 = REAL(getListElement(RsinusoidalForcing,
                                            "t0"))[0];
        sinusoidal_period = REAL(getListElement(RsinusoidalForcing,
                                                "period"))[0];
    } else
    {
        sinusoidal = 0;
    }

    SEXP Ruptake = getListElement(gparms, "vaccine.uptake");
    if (Ruptake != NULL)
    {
        uptake = REAL(Ruptake);
        uptakeTimes = LENGTH(Ruptake) / nAgeGroups;
    } else
    {
        uptakeTimes = 0;
    }

    SEXP Rcampaign = getListElement(gparms, "vaccine.campaign");
    if (Rcampaign != NULL)
    {
        campaign = REAL(Rcampaign);
        campaignTimes = LENGTH(Rcampaign) / nAgeGroups;
    } else
    {
        campaignTimes = 0;
    }

    SEXP RboosterCampaign = getListElement(gparms, "booster.campaign");
    if (RboosterCampaign != NULL)
    {
        boosterCampaign = REAL(RboosterCampaign);
        boosterCampaignTimes = LENGTH(RboosterCampaign) / nAgeGroups;
    } else
    {
        boosterCampaignTimes = 0;
    }

    SEXP Rbeds = getListElement(gparms, "beds");
    if (Rbeds != NULL)
    {
        beds = INTEGER(Rbeds);
        bedTimes = LENGTH(Rbeds) / nAgeGroups;
    } else
    {
        bedTimes = 0;
    }

    started = 0;

    if (verbose)
    {
        if (Rbeta != NULL)
        {
            printf("beta:%f\n", beta);
        }
        if (Rgamma != NULL)
        {
            printf("gamma:%f\n", gamma);
        }
        if (Rdelta != NULL)
        {
            printf("delta:%f\n", delta);
        }
        if (Repsilon != NULL)
        {
            printf("epsilon:%f\n", epsilon);
        }
        if (Rrho != NULL)
        {
            printf("rho:%f\n", rho);
        }
        if (Rasymp != NULL)
        {
            printf("asymp:%f\n", asymp);
        }
        if (Rreservoir != NULL)
        {
            printf("reservoir:%f\n", r);
        }
        if (Rmu != NULL)
        {
            printf("birthrate:%f\n", birthrate);
            printf("deathrate:%f\n", deathrate);
        }
        if (RnE != NULL)
        {
            printf("nE:%d\n", nE);
        }
        if (RnC != NULL)
        {
            printf("nC:%d\n", nC);
        }
        if (RnI != NULL)
        {
            printf("nI:%d\n", nI);
        }

        if (Refficacy != NULL)
        {
            printf("vaccine_efficacy:%f\n", vaccine_efficacy);
        }
        if (Rmaternal_immunity != NULL)
        {
            printf("maternal_immunity:%f\n", maternal_immunity);
        }
        if (Rantibody_development != NULL)
        {
            printf("antibody_development:%f\n", antibody_development);
        }
        if (Rleaky_immunity != NULL)
        {
            printf("leaky_immunity:%f\n", leaky_immunity);
        }
        if (Rboosted_immunity != NULL)
        {
            printf("boosted_immunity:%f\n", boosted_immunity);
        }
        if (Rvaccine_waning != NULL)
        {
            printf("vaccine_waning:%f\n", vaccine_waning);
        }
        if (Rboosted_waning != NULL)
        {
            printf("boosted_waning:%f\n", boosted_waning);
        }
        if (Rhospital_entry != NULL)
        {
            printf("hospital_entry:%f\n", hospital_entry);
        }
        if (Rinterpolate == NULL)
        {
            if (interpolate)
            {
                printf("interpolate: YES\n");
            } else
            {
                printf("interpolate: NO\n");
            }
        }

        if (Rpd != NULL)
        {
            for (unsigned int i = 0; i < nAgeGroups; ++i)
            {
                printf("pd[%u]:%f\n", i, pd[i]);
            }
        }

        for (unsigned int i = 0; i < nAgeGroups; ++i)
        {
            for (unsigned int j = 0; j < nAgeGroups; ++j)
            {
                printf("mixing[%u,%u]:%f\n", i, j, NOR_MIX(i, j));
            }
        }
        if (termtime)
        {
            printf("holiday_start:%f\n", termtime_holiday_start);
            printf("holiday_end:%f\n", termtime_holiday_end);
            for (unsigned int i = 0; i < nAgeGroups; ++i)
            {
                for (unsigned int j = 0; j < nAgeGroups; ++j)
                {
                    printf("holiday_mixing[%u,%u]:%f\n", i, j, HOL_MIX(i, j));
                }
            }
        }
        if (sinusoidal)
        {
            printf("sinusoidal_amplitude:%f\n", sinusoidal_amplitude);
            printf("sinusoidal_t0:%f\n", sinusoidal_t0);
            printf("sinusoidal_period:%f\n", sinusoidal_period);
        }
        for (unsigned int i = 0; i < uptakeTimes; ++i)
        {
            for (unsigned int j = 0; j < nAgeGroups; ++j)
            {
                printf("uptake[%u,%u]:%f\n", i, j,
                       get_fvalues(uptake, i, j,
                                   uptakeTimes, nAgeGroups));
            }
        }
        for (unsigned int i = 0; i < campaignTimes; ++i)
        {
            for (unsigned int j = 0; j < nAgeGroups; ++j)
            {
                printf("campaign[%u,%u]:%f\n", i, j,
                       get_fvalues(campaign, i, j,
                                   campaignTimes, nAgeGroups));
            }
        }
        for (unsigned int i = 0; i < boosterCampaignTimes; ++i)
        {
            for (unsigned int j = 0; j < nAgeGroups; ++j)
            {
                printf("boosterCampaign[%u,%u]:%f\n", i, j,
                       get_fvalues(boosterCampaign, i, j,
                                   boosterCampaignTimes, nAgeGroups));
            }
        }
        for (unsigned int i = 0; i < bedTimes; ++i)
        {
            printf("beds[%u]:%d\n", i,
                   get_ivalue(beds, i, bedTimes));
        }
        for (unsigned int i = 0; i < birthTimes; ++i)
        {
            printf("births[%u]:%d\n", i,
                   get_ivalue(births, i, birthTimes));
        }
        for (unsigned int i = 0; i < deathTimes; ++i)
        {
            printf("deaths[%u]:%d\n", i,
                   get_ivalue(deaths, i, deathTimes));
        }
        for (unsigned int i = 0; i < childMixingTimes; ++i)
        {
            printf("childMixing[%u]:%f\n", i,
                   get_fvalue(childMixing, i, childMixingTimes));
        }
        if (RchildMixingGroup != NULL)
        {
            printf("childMixingGroup:%d\n", childMixingGroup);
        }
        for (unsigned int i = 0; i < (nAgeGroups - 1); ++i)
        {
            printf("ageing[%u]: %f\n", i, ageing[i]);
        }
    }

}

/* derivatives and output variables */
void seir_derivs (int *neq, double *t, double *y, double *ydot,
                  double *yout, int *ip)
{
    if (ip[0] <1) error("nout should be at least 1");

    double change = 0;
    int holidays = 0;
    double time = *t;
    if (!started)
    {
        started = 1;
        start_step = floor(time);
    }

    double integer;
    if (termtime)

    {
        double season = modf(time, &integer);
        if (season >= termtime_holiday_start &&
            season < termtime_holiday_end)

        {
            holidays = 1;
        }
    }

    double Isum[nAgeGroups];
    double Csum[nAgeGroups];
    for (unsigned int i = 0; i < nAgeGroups; ++i)
    {
        Isum[i] = 0;
        for (unsigned int j = 0; j < nI; ++j)

        {
            Isum[i] += I(i, j);
        }
    }
    for (unsigned int i = 0; i < nAgeGroups; ++i)
    {
        Csum[i] = 0;
        for (unsigned int j = 0; j < nC; ++j)

        {
            Csum[i] += C(i, j);
        }
    }

    if (verbose)
    {
        printf("time: %f\n", time);
        if (birthTimes) {
            printf("current births: %d\n", get_ivalue(births, time, birthTimes));
        }
        if (deathTimes) {
            printf("current deaths: %d\n", get_ivalue(deaths, time, deathTimes));
        }
        if (uptakeTimes) {
            printf("current uptake: %f", get_fvalues(uptake, time, 0,
                                                     uptakeTimes, nAgeGroups));
            for (unsigned int i = 1; i < nAgeGroups; ++i)
            {
                printf(",%f", get_fvalues(uptake, time, i,
                                          uptakeTimes, nAgeGroups));
            }
            printf("\n");
        }
        if (campaignTimes) {
            printf("current campaign: %f", get_fvalues(campaign, time, 0,
                                                       campaignTimes, nAgeGroups));
            for (unsigned int i = 1; i < nAgeGroups; ++i)
            {
                printf(",%f", get_fvalues(campaign, time, i,
                                          campaignTimes, nAgeGroups));
            }
            printf("\n");
        }
        if (boosterCampaignTimes) {
            printf("current booster campaign: %f",
                   get_fvalues(boosterCampaign, time, 0,
                               boosterCampaignTimes, nAgeGroups));
            for (unsigned int i = 1; i < nAgeGroups; ++i)
            {
                printf(",%f", get_fvalues(boosterCampaign, time, i,
                                          boosterCampaignTimes, nAgeGroups));
            }
            printf("\n");
        }
    }

    totalN = 0;
    for (unsigned int i = 0; i < nAgeGroups; ++i)
    {
        N(i) = S(i) + Isum[i] + Csum[i] + R(i);
        double Esum = 0;
        for (unsigned int j = 0; j < nE; ++j)
        {
            Esum += E(i, j);
        }
        N(i) += Esum;
        if (maternal_immunity > 0)
        {
            if (i == 0)
            {
                N(i) += B;
            }
            if (verbose)
            {
                printf("B[%u]:%f ", i, B);
            }
        }
        if (verbose)
        {
            printf("S[%u]:%f", i, S(i));
            if (nE > 0)
            {
                printf(" E[%u]:%f", i, Esum);
            }
            if (nC > 0)
            {
                printf(" C[%u]:%f", i, Csum[i]);
            }
            printf(" I[%u]:%f R[%u]:%f",
                   i, Isum[i], i, R(i));
        }
        if (antibody_development > 0)
        {
            N(i) += A(i);
            if (verbose)
            {
                printf(" A[%u]:%f", i, A(i));
            }
        }
        if (uptakeTimes + campaignTimes > 0)
        {
            N(i) += V(i);
            if (verbose)
            {
                printf(" V[%u]:%f", i, V(i));
            }
        }
        if (boosterCampaignTimes > 0)
        {
            N(i) += W(i);
            if (verbose)
            {
                printf(" W[%u]:%f", i, W(i));
            }
        }
        if (verbose)
        {
            printf(" N[%u]:%f\n", i, N(i));
        }
        totalN += N(i);
    }

    // work out death normalisation
    double dn = 0;
    for (unsigned int i = 0; i < nAgeGroups; ++i)
    {
        dn += pd[i] * N(i);
    }

    // births
    double newBirths = 0;

    if (birthTimes > 0)
    {
        newBirths = get_ivalue(births, time, birthTimes);
    } else
    {
        newBirths = birthrate * totalN;
    }

    /* births into susceptibles and vaccination */
    Sdot(0) = (1 - maternal_immunity) *
        (1 - get_fvalues(uptake, time, 0, uptakeTimes, nAgeGroups)) * newBirths;
    if (verbose > 1 && Sdot(0) != 0)
    {
        printf("Sdot[%d] (births): %f\n", 0, Sdot(0));
    }
    /* vaccination in first age group */
    double vaccination =
        get_fvalues(uptake, time, 0, uptakeTimes, nAgeGroups) * newBirths +
        get_fvalues(campaign, time, 0, campaignTimes, nAgeGroups) * newBirths;
    if (antibody_development)
    {
        Adot(0) = vaccination;
        if (verbose > 1 && Adot(0) != 0)
        {
            printf("Adot[%d] (vaccination): %f\n", 0, Adot(0));
        }
    } else if (uptakeTimes > 0)
    {
        Vdot(0) = vaccination;
        if (verbose > 1 && Vdot(0) != 0)
        {
            printf("Vdot[%d] (vaccination): %f\n", 0, Vdot(0));
        }
    }

    /* maternal immune */
    if (maternal_immunity > 0)
    {
        /* births */
        Bdot = maternal_immunity *
            (1 - get_fvalues(uptake, time, 0, uptakeTimes, nAgeGroups)) *
            newBirths;
        if (verbose > 1 && Bdot != 0)
        {
            printf("Bdot (maternal immunity): %f\n", Bdot);
        }
        /* ageing */
        change = - ageing[0] * B;
        Bdot += change;
        if (verbose > 1 && change != 0)
        {
            printf("Bdot (ageing): %f\n", Bdot);
        }
        /* ageing into susceptibles */
        Sdot(1) = (1 - get_fvalues(uptake, time, 1, uptakeTimes, nAgeGroups)) *
            ageing[0] * B;
        if (verbose > 1 && Sdot(1) != 0)
        {
            printf("Sdot[%d] (ageing from maternally immune): %f\n", 1, Bdot);
        }
        /* vaccination into second age group */
        if (uptakeTimes > 0) {
            if (antibody_development)
            {
                Adot(1) =
                    get_fvalues(uptake, time, 1, uptakeTimes, nAgeGroups) *
                    ageing[0] * B;
                if (verbose > 1 && Adot(1) != 0)
                {
                    printf("Adot[%d] (vaccination into second age group): %f\n",
                           1, Adot(1));
                }
            } else
            {
                Vdot(1) =
                    get_fvalues(uptake, time, 1, uptakeTimes, nAgeGroups) *
                    ageing[0] * B;
                if (verbose > 1 && Vdot(1) != 0)
                {
                    printf("Vdot[%d] (vaccination into second age group): %f\n",
                           1, Vdot(1));
                }
            }
        }
    } else
    {
        Sdot(1) = 0;
        if (antibody_development)
        {
            Adot(1) = 0;
        } else if (uptakeTimes > 0) {
            Vdot(1) = 0;
        }
    }

    for (unsigned int i = 0; i < nAgeGroups; ++i)
    {
        if (i > 1)
        {
            Sdot(i) = 0;
        }
        /* ===== susceptibles ===== */
        /* Ageing */
        change = - (i < (nAgeGroups - 1) ? ageing[i] * S(i) : 0);
        Sdot(i) += change;
        if (verbose > 1 && change != 0)
        {
            printf("Sdot[%u] (ageing): %f\n", i, change);
        }
        /* Vaccination */
        change = (1 - get_fvalues(uptake, time, i, uptakeTimes, nAgeGroups)) *
            (i > 0 ? ageing[i - 1] * S(i - 1) : 0);
        Sdot(i) += change;
        if (verbose > 1 && change != 0)
        {
            printf("Sdot[%u] (ageing & vaccination): %f\n", i, change);
        }
        /* Deaths */
        change = - (deathTimes > 0 ? get_ivalue(deaths, time, deathTimes) *
                    pd[i] * S(i) / dn : 0);
        Sdot(i) += change;
        if (verbose > 1 && change != 0)
        {
            printf("Sdot[%u] (deaths): %f\n", i, change);
        }


        /* Force of infection: constant reservoir */
        double foi = r;
        if (verbose)
        {
            printf("foi(%u,r):%f\n", i, r);
        }

        /* Force of infection: mixing */
        for (unsigned int j = 0; j < nAgeGroups; ++j)
        {
            double mult = 1;
            if (i == childMixingGroup && j == childMixingGroup)
            {
                mult = get_fvalue(childMixing, time, childMixingTimes);
            }
            if (holidays == 0)
            {
                foi += mult * beta * NOR_MIX(i, j) * (Isum[j] + Csum[j]) / N(j);
            } else
            {
                foi += mult * beta * HOL_MIX(i, j) * (Isum[j] + Csum[j]) / N(j);
            }
            if (verbose)
            {
                printf("foi(%u,%u):%f\n", i, j, foi);
            }
        }
        /* Sinusoidal forcing */
        if (sinusoidal)
        {
            foi = foi * sinusoidal_amplitude *
                sin(2 * PI * (time - sinusoidal_t0) / sinusoidal_period);
        }

        /* Infection */
        change = - foi * S(i);
        Sdot(i) += change;
        if (verbose > 1 && change != 0)
        {
            printf("Sdot[%u] (infection): %f\n", i, change);
        }

        /* ===== Exposed ===== */
        /* Record change in the first infection class (for incidence)  */
        double changeFirstI = 0;
        if (nE + nC == 0)
        {
            /* change in first infected class (from infection, if
               there is no exposed or class) */
            changeFirstI = (1 - asymp) * foi * S(i);
        } else
        {
            if (nE > 0)
            {
                /* Demography */
                for (unsigned int j = 0; j < nE; ++j)
                {
                    Edot(i, j) = 0;
                    /* Ageing */
                    change =
                        - (i < (nAgeGroups - 1) ? ageing[i] * E(i, j) : 0) +
                        (i > 0 ? ageing[i - 1] * E(i - 1, j) : 0);
                    Edot(i, j) += change;
                    if (verbose > 1 && change != 0)
                    {
                        printf("Edot[%u,%u] (ageing): %f\n", i, j, change);
                    }

                    /* Deaths */
                    change =
                        - (deathTimes > 0 ? get_ivalue(deaths, time, deathTimes) *
                           E(i, j) * pd[i] / dn : 0);
                    Edot(i, j) += change;
                    if (verbose > 1 && change != 0)
                    {
                        printf("Edot[%u,%u] (deaths): %f\n", i, j, change);
                    }
                }
                /* Infection */
                change = foi * S(i);
                Edot(i, 0) += change;
                if (verbose > 1 && change != 0)
                {
                    printf("Edot[%u,%d] (infection): %f\n", i, 0, change);
                }

                /* Moving through exposed classes */
                change = - rho * nE * E(i, 0);
                Edot(i,0) += change;
                if (verbose > 1 && change != 0)
                {
                    printf("Edot[%u,%u] (moving through exposed classes): %f\n",
                           i, 0, change);
                }
                for (unsigned int j = 1; j < nE; ++j)
                {
                    change = rho * nE * E(i, j - 1) - rho * nE * E(i, j);
                    Edot(i, j) += change;
                    if (verbose > 1 && change != 0)
                    {
                        printf("Edot[%u,%u] (moving through exposed classes): %f\n",
                               i, j, change);
                    }
                }
                /* change in first infected class (from incubation, if
                   there is an exposed class) */
                if (nC == 0)
                {
                    changeFirstI = (1 - asymp) * rho * nE * E(i, nE - 1);
                }
            }
            if (nC > 0)
            {
                /* Demography */
                for (unsigned int j = 0; j < nC; ++j)
                {
                    Cdot(i, j) = 0;
                    /* Ageing */
                    change = - (i < (nAgeGroups - 1) ? ageing[i] * C(i, j) : 0) +
                        (i > 0 ? ageing[i - 1] * C(i - 1, j) : 0);
                    Cdot(i, j) += change;
                    if (verbose > 1 && change != 0)
                    {
                        printf("Cdot[%u,%u] (ageing): %f\n", i, j, change);
                    }

                    /* Deaths */
                    change = - (deathTimes > 0 ? get_ivalue(deaths, time, deathTimes) *
                                C(i, j) * pd[i] / dn : 0);
                    Cdot(i, j) += change;
                    if (verbose > 1 && change != 0)
                    {
                        printf("Cdot[%u,%u] (deaths): %f\n", i, j, change);
                    }
                }
                /* finishing incubation */
                change = rho * (1 - asymp) *  nE * E(i, nE - 1);
                Cdot(i, 0) += change;
                if (verbose)
                {
                    printf("Cdot[%u,%d] (finished incubation): %f\n", i, 0, change);
                }

                /* Moving through community classes */
                change = - epsilon * nC * C(i, 0);
                Cdot(i,0) += change;
                if (verbose > 1 && change != 0)
                {
                    printf("Cdot[%u,%u] (moving through community classes): %f\n",
                           i, 0, change);
                }
                for (unsigned int j = 1; j < nC; ++j)
                {
                    change = epsilon * nC * C(i, j - 1) - epsilon * nC * C(i, j);
                    Cdot(i, j) += change;
                    if (verbose > 1 && change != 0)
                    {
                        printf("Cdot[%u,%u] (moving through community classes): %f\n",
                               i, j, change);
                    }
                }
                /* change in first infected class (from community, if
                   there is an community class) */
                changeFirstI = epsilon * nC * C(i, nC - 1);
            }
        }

        /* ===== vaccination ===== */
        if (uptakeTimes + campaignTimes > 0)
        {
            if (i > 1 | antibody_development > 0)
            {
                Vdot(i) = 0;
            }
        }
        if (antibody_development > 0)
        {
            if (i > 1)
            {
                Adot(i) = 0;
            }

            /* Ageing */
            change = - (i < (nAgeGroups - 1) ? ageing[i] * A(i) : 0) +
                (i > 0 ? ageing[i - 1] * A(i - 1) : 0);
            Adot(i) += change;
            if (verbose > 1 && change != 0)
            {
                printf("Adot[%u] (ageing): %f\n", i, change);
            }
            /* Vaccination */
            change = get_fvalues(uptake, time, i, uptakeTimes, nAgeGroups) *
                (i > 0 ? ageing[i - 1] * S(i - 1) : 0);
            Adot(i) += change;
            if (verbose > 1 && change != 0)
            {
                printf("Adot[%u] (vaccination): %f\n", i, change);
            }
            /* Campaign vaccination */
            change = fmax(fmin(S(i) - Sdot(i),
                               get_fvalues(campaign, time, i, campaignTimes,
                                           nAgeGroups)),
                          0);

            Adot(i) += change;
            Sdot(i) -= change;
            if (verbose > 1 && change != 0)
            {
                printf("Adot[%u] (campaign vaccination): %f\n", i, change);
                printf("Sdot[%u] (campaign vaccination): %f\n", i, - change);
            }
            /* Infection */
            change = - foi * A(i);
            Adot(i) += change;
            if (verbose)
            {
                printf("Adot[%u] (infection): %f\n", i, change);
            }
            if (nE > 0)
            {
                Edot(i, 0) -= change;
                if (verbose > 1 && change != 0)
                {
                    printf("Edot[%u,%u] (infection from building antibodies): %f\n", i, 0,
                           - change);
                }
            } else
            {
                changeFirstI -= change;
                if (verbose > 1 && change != 0)
                {
                    printf("Idot[%u,%u] (infection from building antibodies): %f\n", i, 0,
                           - change);
                }
            }
            /* Maturing antibodies */
            change = - antibody_development * A(i);
            Adot(i) += change;
            if (verbose > 1 && change != 0)
            {
                printf("Adot[%u] (maturing antibodies): %f\n", i, change);
            }
            /* death */
            change = - (deathTimes > 0 ? get_ivalue(deaths, time, deathTimes) *
                        A(i) * pd[i] / dn : 0);
            Adot(i) += change;
            if (verbose > 1 && change != 0)
            {
                printf("Adot[%u] (death): %f\n", i, change);
            }
        }
        if (uptakeTimes + campaignTimes > 0)
        {
            /* Ageing */
            change = - (i < (nAgeGroups - 1) ? ageing[i] * V(i) : 0) +
                (i > 0 ? ageing[i - 1] * V(i - 1) : 0);
            Vdot(i) += change;
            if (verbose > 1 && change != 0)
            {
                printf("Vdot[%u] (ageing): %f\n", i, change);
            }
            /* Vaccination */
            if (antibody_development == 0)
            {
                change = get_fvalues(uptake, time, i, uptakeTimes, nAgeGroups) *
                    (i > 0 ? ageing[i - 1] * S(i - 1) : 0);
                Vdot(i) += change;
                if (verbose > 1 && change != 0)
                {
                    printf("Vdot[%u] (vaccination): %f\n", i, change);
                }
                /* Campaign vaccination */
                change = fmax(fmin(S(i) - Sdot(i),
                                   get_fvalues(campaign, time, i, campaignTimes,
                                               nAgeGroups)),
                              0);

                Vdot(i) += change;
                Sdot(i) -= change;
                if (verbose > 1 && change != 0)
                {
                    printf("Vdot[%u] (campaign vaccination): %f\n", i, change);
                    printf("Sdot[%u] (campaign vaccination): %f\n", i, - change);
                }
            } else
            {
                /* Maturing antibodies */
                change = antibody_development * A(i);
                Vdot(i) += change;
                if (verbose > 1 && change != 0)
                {
                    printf("Vdot[%u] (maturing antibodies): %f\n", i, change);
                }
            }
            /* Infection */
            change = - foi * (1 - leaky_immunity) * V(i);
            Vdot(i) += change;
            if (verbose > 1 && change != 0)
            {
                printf("Vdot[%u] (infection): %f\n", i, change);
            }
            if (nE > 0)
            {
                Edot(i, 0) -= change;
                if (verbose > 1 && change != 0)
                {
                    printf("Edot[%u,%u] (infection from vaccinated): %f\n", i, 0,
                           - change);
                }
            } else
            {
                changeFirstI -= change;
                if (verbose > 1 && change != 0)
                {
                    printf("Idot[%u,%u] (infection from vaccinated): %f\n", i, 0,
                           - change);
                }
            }
            /* Waning */
            change = - vaccine_waning * V(i);
            Vdot(i) += change;
            Sdot(i) -= change;

            if (verbose > 1 && change != 0)
            {
                printf("Vdot[%u] (waning immunity): %f\n", i, change);
                printf("Sdot[%u] (waning immunity): %f\n", i, - change);
            }
            /* Death */
            change = - (deathTimes > 0 ? get_ivalue(deaths, time, deathTimes) *
                        V(i) * pd[i] / dn : 0);
            Vdot(i) -= change;
            if (verbose > 1 && change != 0)
            {
                printf("Adot[%u] (death): %f\n", i, change);
            }
        }
        if (boosterCampaignTimes > 0)
        {
            Wdot(i) = 0;
            change = - fmax(fmin(V(i) - Vdot(i),
                                 get_fvalues(boosterCampaign, time, i,
                                             boosterCampaignTimes, nAgeGroups)),
                            0);
            Vdot(i) -= change;
            Wdot(i) += change;
            if (verbose > 1 && change != 0)
            {
                printf("Vdot[%u] (campaign boosting): %f\n", i, change);
                printf("Wdot[%u] (campaign boosting): %f\n", i, change);
            }
            /* Infection */
            change = - foi * (1 - boosted_immunity) * W(i);
            Wdot(i) += change;
            if (verbose > 1 && change != 0)
            {
                printf("Wdot[%u] (infection): %f\n", i, change);
            }
            if (nE > 0)
            {
                Edot(i, 0) -= change;
                if (verbose > 1 && change != 0)
                {
                    printf("Edot[%u,%u] (infection from vaccinated): %f\n", i, 0,
                           - change);
                }
            } else
            {
                changeFirstI -= change;
                if (verbose)
                {
                    printf("Idot[%u,%u] (infection from vaccinated): %f\n", i, 0,
                           - change);
                }
            }
            /* Waning */
            change = - boosted_waning * W(i);
            Wdot(i) += change;
            Sdot(i) -= change;

            if (verbose > 1 && change != 0)
            {
                printf("Wdot[%u] (waning immunity): %f\n", i, change);
                printf("Sdot[%u] (waning immunity): %f\n", i, change);
            }
            /* Death */
            change = - (deathTimes > 0 ? get_ivalue(deaths, time, deathTimes) *
                        W(i) * pd[i] / dn : 0);
            Wdot(i) += change;
            if (verbose > 1 && change != 0)
            {
                printf("Wdot[%u] (death): %f\n", i, change);
            }
        }

        /* ===== infected ===== */
        /* Demography */
        for (unsigned int j = 0; j < nI; ++j)
        {
            Idot(i, j) = 0;
            /* Ageing */
            change = - (i < (nAgeGroups - 1) ? ageing[i] * I(i, j) : 0) +
                (i > 0 ? ageing[i - 1] * I(i - 1, j) : 0);

            Idot(i, j) += change;
            if (verbose > 1 && change != 0)
            {
                printf("Idot[%u,%u] (ageing): %f\n", i, j, change);
            }
            /* Deaths */
            change = - (deathTimes > 0 ? get_ivalue(deaths, time, deathTimes) *
                        I(i, j) * pd[i] / dn : 0);
            Idot(i, j) += change;

            if (verbose > 1 && change != 0)
            {
                printf("Idot[%u,%u] (deaths): %f\n", i, j, change);
            }
        }
        /* Becoming infected & recovery */
        change = changeFirstI;
        Idot(i, 0) += change;
        if (verbose > 1 && change != 0)
        {
            printf("Idot[%u,%d] (infection): %f\n", i, 0, change);
        }
        change = - gamma * nI * I(i, 0);
        Idot(i, 0) += change;
        if (verbose > 1 && change != 0)
        {
            printf("Idot[%u,%u] (moving through infected classes): %f\n",
                   i, 0, change);
        }
        /* Moving through infected classes */
        for (unsigned int j = 1; j < nI; ++j)
        {
            change = gamma * nI * I(i, j - 1) - gamma * nI * I(i, j);
            Idot(i, j) += change;
            if (verbose > 1 && change != 0)
            {
                printf("Idot[%u,%u] (moving through infected classes): %f\n",
                       i, j, change);
            }
        }

        /* ===== recovered ===== */
        Rdot(i) = 0;
        /* Ageing */
        change = - (i < (nAgeGroups - 1) ? ageing[i] * R(i) : 0) +
            (i > 0 ? ageing[i - 1] * R(i - 1) : 0);
        Rdot(i) += change;
        if (verbose > 1 && change != 0)
        {
            printf("Rdot[%u] (ageing): %f\n", i, change);
        }
        /* asymptomatic infection */
        if (nE == 0)
        {
            change = asymp * foi * S(i);
        } else
        {
            change = rho * asymp * nE * E(i, nE - 1);
        }
        Rdot(i) += change;
        if (verbose > 1 && change != 0)
        {
            printf("Rdot[%u] (asymptomatic infection): %f\n", i, change);
        }
        /* Deaths */
        change = - (deathTimes > 0 ? get_ivalue(deaths, time, deathTimes) *
                    R(i) * pd[i] / dn : 0);
        Rdot(i) += change;
        if (verbose > 1 && change != 0)
        {
            printf("Rdot[%u] (deaths): %f\n", i, change);
        }
        /* Recovery */
        change = gamma * nI * I(i, nI - 1);
        Rdot(i) += change;
        if (verbose > 1 && change != 0)
        {
            printf("Rdot[%u] (recovery): %f\n", i, change);
        }
        /* Incidence */
        Zdot(i) = changeFirstI;
        if (verbose > 1 && changeFirstI > 0)
        {
            printf("Zdot[%u] (incidence): %f\n", i, changeFirstI);
        }
    }

    /* deaths */
    if (deathTimes == 0)
    {
        for (unsigned int i = 0; i < nAgeGroups; ++i)
        {
            if (i == 0)
            {
                change = - deathrate * B;
                Bdot += change;
                if (verbose > 1 && change != 0)
                {
                    printf("Bdot (deaths): %f\n", change);
                }
            }
            change = - deathrate * S(i);
            Sdot(i) += change;
            if (verbose > 1 && change != 0)
            {
                printf("Sdot[%u] (deaths): %f\n", i, change);
            }
            for (unsigned int j = 0; j < nE; ++j)
            {
                change = - deathrate * E(nAgeGroups - 1, j);
                Edot(i, j) += change;
                if (verbose > 1 && change != 0)
                {
                    printf("Edot[%u,%u] (deaths): %f\n", i, j, change);
                }
            }
            for (unsigned int j = 0; j < nI; ++j)
            {
                change = - deathrate * C(nAgeGroups - 1, j);
                Cdot(i, j) += change;
                if (verbose > 1 && change != 0)
                {
                    printf("Cdot[%u,%u] (deaths): %f\n", i, j, change);
                }
            }
            for (unsigned int j = 0; j < nI; ++j)
            {
                change = - deathrate * I(nAgeGroups - 1, j);
                Idot(i, j) += change;
                if (verbose > 1 && change != 0)
                {
                    printf("Idot[%u,%u] (deaths): %f\n", i, j, change);
                }
            }
            if (uptakeTimes + campaignTimes > 0)
            {
                change = - deathrate * V(i);
                Vdot(i) -= change;
                if (verbose > 1 && change != 0)
                {
                    printf("Vdot[%u] (deaths): %f\n", i, change);
                }
            }
            if (boosterCampaignTimes > 0)
            {
                change = - deathrate * W(i);
                Wdot(i) -= change;
                if (verbose > 1 && change != 0)
                {
                    printf("Wdot[%u] (deaths): %f\n", i, change);
                }
            }
            if (antibody_development > 0)
            {
                change = - deathrate * A(i);
                Adot(i) -= change;
                if (verbose > 1 && change != 0)
                {
                    printf("Adot[%u] (deaths): %f\n", i, change);
                }
            }
            Rdot(i) -= deathrate * R(i);
            if (verbose)
            {
                printf("Rdot[%u] (deaths): %f\n", i, change);
            }
        }

    }

    if (verbose)
    {
        for (unsigned int i = 0; i < nAgeGroups; ++i)
        {
            double Ndot = Sdot(i);
            printf("Sdot[%u]:%f", i, Sdot(i));
            if (nE > 0)
            {
                double EdotSum = 0;
                for (unsigned int j = 0; j < nE; ++j)
                {
                    EdotSum += Edot(i, j);
                }
                Ndot += EdotSum;
                printf(" Edot[%u]:%f", i, EdotSum);
            }
            if (nC > 0)
            {
                double CdotSum = 0;
                for (unsigned int j = 0; j < nE; ++j)
                {
                    CdotSum += Cdot(i, j);
                }
                Ndot += CdotSum;
                printf(" Cdot[%u]:%f", i, CdotSum);
            }
            double IdotSum = 0;
            for (unsigned int j = 0; j < nE; ++j)
            {
                IdotSum += Idot(i, j);
            }
            Ndot += IdotSum + Rdot(i);
            printf(" Idot[%u]:%f Rdot[%u]:%f",
                   i, IdotSum, i, Rdot(i));
            if (antibody_development > 0)
            {
                Ndot += Adot(i);
                printf(" Adot[%u]:%f", i, Adot(i));
            }
            if (uptakeTimes + campaignTimes > 0)
            {
                Ndot += Vdot(i);
                printf(" Vdot[%u]:%f", i, Vdot(i));
            }
            if (boosterCampaignTimes > 0)
            {
                Ndot += Wdot(i);
                printf(" Wdot[%u]:%f", i, Wdot(i));
            }
            printf(" Ndot[%u]:%f\n", i, Ndot);
        }
    }

    /* fix numerical instabilities */
    for (unsigned int i = 0; i < *neq; ++i)
    {
        if (y[i] < 1 && y[i] + ydot[i] < 0)
        {
            ydot[i] = -y[i];
        }
    }
}
