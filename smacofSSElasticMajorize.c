#include "smacofSSElastic.h"

void smacofSSElasticMajorize(int* nobj, int* ndim, int* ndat, double* snew,
                             int* iind, int* jind, int* weighted, double* wght,
                             double* vinv, double* edis, double* dhat,
                             double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    double* xtmp = xmalloc(Nobj * Ndim * sizeof(double));
    for (int k = 0; k < Nobj * Ndim; k++) {
        xtmp[k] = 0.0;
    }
    for (int k = 0; k < Ndat; k++) {
        if (edis[k] == 0.0) {
            continue;
        }
        int is = iind[k], js = jind[k];
        double elem = wght[k] * dhat[k] / edis[k];
        for (int s = 0; s < Ndim; s++) {
            double add = elem * (xold[is] - xold[js]);
            xtmp[is] += add;
            xtmp[js] -= add;
            is += Nobj;
            js += Nobj;
        }
    }
    for (int k = 0; k < Nobj * Ndim; k++) {
        if (weighted) {
            xnew[k] = 0.0;
        } else {
            xnew[k] = xtmp[k];
        }
    }
    if (weighted) {
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
    } else {
        for (int i = 0; i < Nobj; i++) {
            int is = i;
            for (int s = 0; s < Ndim; s++) {
                xnew[is] /= (double)Nobj;
                is += Nobj;
            }
        }
    }
    for (int k = 0; k < Ndat; k++) {
        int is = iind[k], js = jind[k];
        double sum = 0.0;
        for (int s = 0; s < Ndim; s++) {
            sum += SQUARE(xnew[is] - xnew[js]);
            is += Nobj;
            js += Nobj;
        }
        edis[k] = sqrt(sum);
    }
    xfree(xtmp);
    return;
}