#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SQUARE(x) ((x) * (x))

void primaryApproach(int *, int *, double *, double *, double *, int *, int *);
void secondaryApproach(int *, int *, double *, double *);
void tertiaryApproach(int *, int *, double *, double *);
void tieBlockAverages(int *, int *, int *, double *, double *, double *,
                      double *, int *, double *);
void monotone(int *, double *, double *);
int myComp(const void *, const void *);
void mySort(double *, double *, double *, int *, int *, int *);

struct quintuple {
    double value;
    double weight;
    double dist;
    int index;
    int jndex;
};

void smacofElasticEngine(int *nobj, int *ndim, int *ndat, int *itel, int *ties,
                      int *itmax, int *digits, int *width, bool *verbose, double *wsum,
                      double *sold, double *snew, double *eps, int *iind,
                      int *jind, int *blks, double *edis, double *dhat,
                      double *wght, double *vinv, double *xold, double *xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    while (true) {
        double *xtmp = (double *)calloc(Nobj * Ndim, sizeof(double));
        for (int k = 0; k < Ndat; k++) {
            if (edis[k] == 0.0) {
                continue;
            }
            int is = iind[k] - 1, js = jind[k] - 1;
            double elem = wght[k] * dhat[k] / edis[k];
            for (int s = 0; s < Ndim; s++) {
                double add = elem * (xold[is] - xold[js]);
                xtmp[is] += add;
                xtmp[js] -= add;
                is += Nobj;
                js += Nobj;
            }
        }
        int k = 0;
        for (int j = 0; j < Nobj - 1; j++) {
            for (int i = j + 1; i < Nobj; i++) {
                double elem = vinv[k];
                int is = i, js = j;
                for (int s = 0; s < Ndim; s++) {
                    double add = elem * (xtmp[is] - xtmp[js]);
                    xnew[is] += add;
                    xnew[js] -= add;
                    is += Nobj;
                    js += Nobj;
                }
                k++;
            }
        }
        free(xtmp);
        for (int k = 0; k < Ndat; k++) {
            int is = iind[k] - 1, js = jind[k] - 1;
            double sum = 0.0;
            for (int s = 0; s < Ndim; s++) {
                sum += SQUARE(xnew[is] - xnew[js]);
                is += Nobj;
                js += Nobj;
            }
            edis[k] = sqrt(sum);
        }
        double smid = 0.0;
        for (int k = 0; k < Ndat; k++) {
            smid += wght[k] * SQUARE(dhat[k] - edis[k]);
        }
        smid /= *wsum;
        double *wwrk = (double *)calloc(Ndat, sizeof(double));
        for (int k = 0; k < Ndat, k++)
        {
            dhat[k] = -1.0 / edis[k];
            wwrk[k] = wght[k] * SQUARE(edis[k])
        }
        if (*ties == 1) {
            (void)primaryApproach(ndat, blks, dhat, wwrk, edis, iind, jind);
        }
        if (*ties == 2) {
            (void)secondaryApproach(ndat, blks, dhat, wwrk);
        }
        if (*ties == 3) {
            (void)tertiaryApproach(ndat, blks, dhat, wwrk);
        }
        double ssq = 0.0;
        for (int k = 0; k < Ndat; k++) {
            dhat[k] = 1 / dhat[k];
            ssq += wght[k] * SQUARE(dhat[k]);
        }
        *snew = 0.0;
        for (int k = 0; k < Ndat; k++) {
            dhat[k] *= sqrt(*wsum / ssq);
            *snew += wght[k] * SQUARE(dhat[k] - edis[k]);
        }
        *snew /= *wsum;
        if (*verbose) {
            printf("itel %4d sold %*.*f smid %*.*f snew %*.*f\n", *itel, *width,
                   *digits, *sold, *width, *digits, smid, *width, *digits,
                   *snew);
        }
        if ((*itel == *itmax) || ((*sold - *snew) < *eps)) {
            break;
        }
        for (int k = 0; k < Nobj * Ndim; k++) {
            xold[k] = xnew[k];
            xnew[k] = 0.0;
        }
        *sold = *snew;
        *itel += 1;
        free(wwrk);
    }
}

void primaryApproach(int *ndat, int *blks, double *x, double *w, double *d,
                     int *iind, int *jind) {
    int Ndat = *ndat;
    for (int i = 0; i < Ndat; i++) {
        int blksize = blks[i];
        if (blksize > 0) {
            double *extracx = (double *)calloc(blksize, sizeof(double));
            double *extracw = (double *)calloc(blksize, sizeof(double));
            double *extracd = (double *)calloc(blksize, sizeof(double));
            int *extraci = (int *)calloc(blksize, sizeof(int));
            int *extracj = (int *)calloc(blksize, sizeof(int));
            for (int j = 0; j < blksize; j++) {
                extracx[j] = x[i + j];
                extracw[j] = w[i + j];
                extracd[j] = d[i + j];
                extraci[j] = iind[i + j];
                extracj[j] = jind[i + j];
            }
            (void)mySort(extracx, extracw, extracd, extraci, extracj, &blksize);
            for (int j = 0; j < blksize; j++) {
                x[i + j] = extracx[j];
                w[i + j] = extracw[j];
                d[i + j] = extracd[j];
                iind[i + j] = extraci[j];
                jind[i + j] = extracj[j];
            }
            free(extracx);
            free(extracw);
            free(extracd);
            free(extraci);
            free(extracj);
        }
    }
    double *ww = (double *)calloc(Ndat, sizeof(double));
    for (int i = 0; i < Ndat; i++) {
        ww[i] = w[i];
    }
    (void)monotone(ndat, x, ww);
    free(ww);
    return;
}

void secondaryApproach(int *ndat, int *blks, double *x, double *w) {
    int nblk = 0, Ndat = *ndat;
    for (int k = 0; k < Ndat; k++) {
        if (blks[k] > 0) {
            nblk++;
        }
    }
    double *xsum = (double *)calloc((size_t)nblk, sizeof(double));
    double *wsum = (double *)calloc((size_t)nblk, sizeof(double));
    double *xave = (double *)calloc((size_t)nblk, sizeof(double));
    int *csum = (int *)calloc((size_t)nblk, sizeof(int));
    (void)tieBlockAverages(ndat, &nblk, blks, x, w, xsum, wsum, csum, xave);
    (void)monotone(&nblk, xave, wsum);
    int l = 1;
    for (int k = 0; k < nblk; k++) {
        for (int i = l; i <= l + csum[k] - 1; i++) {
            x[i - 1] = xave[k];
        }
        l += csum[k];
    }
    free(xsum);
    free(wsum);
    free(csum);
    free(xave);
    return;
}

void tertiaryApproach(int *ndat, int *blks, double *x, double *w) {
    int nblk = 0, Ndat = *ndat;
    for (int k = 0; k < Ndat; k++) {
        if (blks[k] > 0) {
            nblk++;
        }
    }
    double *xsum = (double *)calloc((size_t)nblk, sizeof(double));
    double *wsum = (double *)calloc((size_t)nblk, sizeof(double));
    double *xave = (double *)calloc((size_t)nblk, sizeof(double));
    double *yave = (double *)calloc((size_t)nblk, sizeof(double));
    int *csum = (int *)calloc((size_t)nblk, sizeof(int));
    (void)tieBlockAverages(ndat, &nblk, blks, x, w, xsum, wsum, csum, xave);
    for (int k = 0; k < nblk; k++) {
        yave[k] = xave[k];
    }
    (void)monotone(&nblk, xave, wsum);
    int l = 1;
    for (int k = 0; k < nblk; k++) {
        for (int i = l; i <= l + csum[k] - 1; i++) {
            x[i - 1] = xave[k] + (x[i - 1] - yave[k]);
        }
        l += csum[k];
    }
    free(xsum);
    free(wsum);
    free(xave);
    free(yave);
    free(csum);
    return;
}

void tieBlockAverages(int *ndat, int *nblk, int *blks, double *x, double *w,
                      double *xsum, double *wsum, int *csum, double *xave) {
    int iblk = 0, Ndat = *ndat, Nblk = *nblk;
    for (int k = 0; k < Ndat; k++) {
        if (blks[k] > 0) {
            double sum1 = 0.0, sum2 = 0.0;
            for (int l = k; l < k + blks[k]; l++) {
                sum1 += w[l] * x[l];
                sum2 += w[l];
            }
            xsum[iblk] = sum1;
            wsum[iblk] = sum2;
            csum[iblk] = blks[k];
            iblk++;
        }
    }
    for (int i = 0; i < Nblk; i++) {
        xave[i] = xsum[i] / wsum[i];
    }
}

int myComp(const void *px, const void *py) {
    double x = ((struct quintuple *)px)->value;
    double y = ((struct quintuple *)py)->value;
    return (int)copysign(1.0, x - y);
}

void mySort(double *x, double *w, double *d, int *iind, int *jind, int *n) {
    int nn = *n;
    struct quintuple *xi = (struct quintuple *)calloc(
        (size_t)nn, (size_t)sizeof(struct quintuple));
    for (int i = 0; i < nn; i++) {
        xi[i].value = x[i];
        xi[i].weight = w[i];
        xi[i].dist = d[i];
        xi[i].index = iind[i];
        xi[i].jndex = jind[i];
    }
    (void)qsort(xi, (size_t)nn, (size_t)sizeof(struct quintuple), myComp);
    for (int i = 0; i < nn; i++) {
        x[i] = xi[i].value;
        w[i] = xi[i].weight;
        d[i] = xi[i].dist;
        iind[i] = xi[i].index;
        jind[i] = xi[i].jndex;
    }
    free(xi);
}

void monotone(int *n, double *x, double *w)
// Function monotone(),
// performs simple linear ordered monotone regression
// Copyright (C) 2020 Frank M.T.A. Busing (e-mail: busing at fsw dot leidenuniv
// dot nl) This function is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation. This program is distributed in the hope that it
// will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
// Public License for more details. You should have received a copy of the GNU
// General Public License along with this function. If not, see
// <https://www.gnu.org/licenses/>.
{
    double *rx = &x[-1];
    double *rw = &w[-1];
    size_t *idx = (size_t *)calloc(*n + 1, sizeof(size_t));
    idx[0] = 0;
    idx[1] = 1;
    size_t b = 1;
    double xbm1 = rx[b];
    double wbm1 = rw[b];
    for (size_t i = 2; i <= *n; i++) {
        b++;
        double xb = rx[i];
        double wb = rw[i];
        if (xbm1 > xb) {
            b--;
            double sb = wbm1 * xbm1 + wb * xb;
            wb += wbm1;
            xb = sb / wb;
            while (i < *n && xb >= rx[i + 1]) {
                i++;
                sb += rw[i] * rx[i];
                wb += rw[i];
                xb = sb / wb;
            }
            while (b > 1 && rx[b - 1] > xb) {
                b--;
                sb += rw[b] * rx[b];
                wb += rw[b];
                xb = sb / wb;
            }
        }
        rx[b] = xbm1 = xb;
        rw[b] = wbm1 = wb;
        idx[b] = i;
    }
    size_t from = *n;
    for (size_t k = b; k > 0; k--) {
        const size_t to = idx[k - 1] + 1;
        const double xk = rx[k];
        for (size_t i = from; i >= to; i--) rx[i] = xk;
        from = to - 1;
    }
    free(idx);
}  // monotone
