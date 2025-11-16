#include "smacofSSElastic.h"

void smacofSSElasticEngine(int* nobj, int* ndim, int* ndat, int* itel,
                           int* ties, int* itmax, int* digits, int* width,
                           int* verbose, int* ordinal, int* weighted,
                           double* sold, double* snew, double* eps, int* iind,
                           int* jind, int* blks, double* wght, double* edis,
                           double* dhat, double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    double* wadj = xmalloc(Ndat * sizeof(double));
    double* vinv = xmalloc(Nobj * (Nobj - 1) * sizeof(double) / 2);
    while (true) {
        for (int k = 0; k < Ndat; k++) {
            wadj[k] = wght[k] / SQUARE(dhat[k]);
        }
        (void)smacofMPInverseV(nobj, ndat, iind, jind, wadj, vinv);
        (void)smacofSSElasticMajorize(nobj, ndim, ndat, snew, iind, jind,
                                      weighted, wadj, vinv, edis, dhat, xold,
                                      xnew);
        double smid = smacofSSElasticLoss(ndat, edis, dhat, wght);
        if (*ordinal) {
            for (int k = 0; k < Ndat; k++) {
                wadj[k] = wght[k] * SQUARE(edis[k]);
            }
            for (int k = 0; k < Ndat; k++) {
                dhat[k] = -1 / edis[k];
            }
            (void)smacofSSElasticMonotone(ndat, ties, snew, iind, jind, blks,
                                          edis, dhat, wadj);
            for (int k = 0; k < Ndat; k++) {
                dhat[k] = -1 / dhat[k];
            }
            *snew = smacofSSElasticLoss(ndat, edis, dhat, wght);
        } else {
            *snew = smid;
        }
        if (*verbose) {
            if (*ordinal) {
                printf("itel %4d sold %*.*f smid %*.*f snew %*.*f\n", *itel,
                       *width, *digits, *sold, *width, *digits, smid, *width,
                       *digits, *snew);
            } else {
                printf("itel %4d sold %*.*f snew %*.*f\n", *itel, *width,
                       *digits, *sold, *width, *digits, *snew);
            }
        }
        if ((*itel == *itmax) || ((*sold - *snew) < *eps)) {
            break;
        }
        for (int k = 0; k < Nobj * Ndim; k++) {
            xold[k] = xnew[k];
        }
        *sold = *snew;
        *itel += 1;
    }
    xfree(vinv);
    xfree(wadj);
    return;
}

double smacofSSElasticLoss(int* ndat, double* edis, double* dhat,
                           double* wght) {
    double loss = 0.0;
    for (int k = 0; k < *ndat; k++) {
        loss += wght[k] * SQUARE(1 - edis[k] / dhat[k]);
    }
    return loss;
}