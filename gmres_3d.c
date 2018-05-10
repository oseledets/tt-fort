/**
 *  \file gmres_3d.c
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// extern void tt_adapt_als_mp_djac_apply_(char *ptype, long *rx1, long *n, long *rx2, double *jacs, double *x, double *y, double *work1);
// extern void tt_adapt_als_mp_dbfun3_(long *rx1, long *m, long *rx2, long *ry1, long *n, long *ry2, long *ra1, long *ra2, double *phi1, double *A, double *phi2, double *x, double *y, double *res1, double *res2);
extern void djac_apply(char *ptype, long *rx1, long *n, long *rx2, double *jacs, double *x, double *y, double *work1);
extern void dbfun3(long *rx1, long *m, long *rx2, long *ry1, long *n, long *ry2, long *ra1, long *ra2, double *phi1, double *A, double *phi2, double *x, double *y, double *res1, double *res2);

extern double dnrm2_(long *n, double *x, long *step);
extern double ddot_(long *n, double *x, long *xstep, double *y, long *ystep);

// #define djac_apply tt_adapt_als_mp_djac_apply_
// #define dbfun3 tt_adapt_als_mp_dbfun3_

char gmres_3d_printf_buf_[128];


void dgmresr_hh(double *Phi1, double *A, double *Phi2, double *rhs, long rx1, long n, long rx2, long ra1, long ra2, long nrestart, double tol, long niters, char ptype, double *jacs, double *sol, char verb)
// right preconditioned - for residual tolerance
// This one with Householder tranforms
{
    long i,j,it, sz;
    double *U, *w;
    double *R, *J, *tau;
    double *res1, *res2;
    double nrmr, curres, nrmrhs;
    double dalpha, dbeta;
    char last_iter = 0;
//     double done=1.0;
//     double dzero=0.0;
    long ione = 1;
    char trans = 'N';
    char uplo = 'U';

//     printf("sz: %dx%dx%d, ra:%dx%d, its: %dx%d, prec: %c, verb:%d, tol:%g, long: %d\n",rx1,n,rx2,ra1,ra2,nrestart,niters, ptype, verb, tol, sizeof(long));

    sz = rx1*n*ra2*rx2;
    if (rx1*ra1*n*rx2>sz) sz = rx1*ra1*n*rx2;
    res1 = (double *)malloc(sizeof(double)*sz);
    res2 = (double *)malloc(sizeof(double)*sz);

    sz = rx1*n*rx2;

//     printf("rhs:\n");
//     for (i=0; i<sz; i++) printf("%g\n", rhs[i]);
//     printf("phi1:\n");
//     for (i=0; i<rx1*ra1*rx1; i++) printf("%g\n", Phi1[i]);


    U = (double *)malloc(sizeof(double)*sz*(nrestart+1));
    w = (double *)malloc(sizeof(double)*sz*2);
//     w2 = (double *)malloc(sizeof(double)*sz);
    tau = (double *)malloc(sizeof(double)*(nrestart+1));
    J = (double *)malloc(sizeof(double)*2*(nrestart));
//     Q = (double *)malloc(sizeof(double)*nrestart*(nrestart+1));
    R = (double *)malloc(sizeof(double)*(nrestart)*(nrestart));

    for (j=0; j<sz; j++) sol[j]=0.0;

    for (it=0; it<niters; it++) {
        // r0
        if ((ptype!='n')&&(ptype!='N')) {
            djac_apply(&ptype, &rx1, &n, &rx2, jacs, sol, w, res1);
            dbfun3(&rx1,&n,&rx2, &rx1,&n,&rx2, &ra1, &ra2, Phi1,A,Phi2, w, w, res1, res2);
        }
        else {
            dbfun3(&rx1,&n,&rx2, &rx1,&n,&rx2, &ra1, &ra2, Phi1,A,Phi2, sol, w, res1, res2);
//             dbfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,sol,w);
        }
        dbeta = -1.0;
        daxpy_(&sz,&dbeta,rhs,&ione,w,&ione);
        dscal_(&sz,&dbeta,w,&ione);
//         for (j=0; j<sz; j++) mexPrintf("%3.5e\n", w[j]);
        nrmr = dnrm2_(&sz,w,&ione);
        if (verb>1) printf("restart %d, res: %3.5e\n", it, nrmr);
//         if (verb>1) sprintf(&gmres_3d_printf_buf_[0], "restart %d, res: %3.5e\n", it, nrmr);
        if (it==0) nrmrhs = nrmr;
        // initial HHT
        dbeta = nrmr;
        if (w[0]<0.0) dbeta = -dbeta;
        w[0]+=dbeta;
        tau[0]=-dbeta;
        nrmr = dnrm2_(&sz,w,&ione);
        dbeta = 1.0/nrmr;
        dscal_(&sz,&dbeta,w,&ione);
        dcopy_(&sz,w,&ione,&U[0],&ione);

        for (j=0; j<nrestart; j++) {
            // HHT on last U
            dcopy_(&sz,&U[sz*j],&ione,w,&ione);
            dbeta = -2.0*U[j+sz*j];
            dscal_(&sz,&dbeta,w,&ione);
            w[j]+=1.0;
            for (i=j-1; i>=0; i--) {
                dbeta = -2.0*ddot_(&sz,&U[i*sz],&ione,w,&ione);
                daxpy_(&sz,&dbeta,&U[i*sz],&ione,w,&ione);
            }
            dbeta = dnrm2_(&sz,w,&ione);
            dbeta = 1.0/dbeta;
            dscal_(&sz,&dbeta,w,&ione); // w=w/norm(w);

            // precvec, matvec
            if ((ptype!='n')&&(ptype!='N')) {
                djac_apply(&ptype, &rx1, &n, &rx2, jacs, w, w, res1);
                dbfun3(&rx1,&n,&rx2, &rx1,&n,&rx2, &ra1, &ra2, Phi1,A,Phi2, w, w, res1, res2);
//                 dcopy_(&sz,w,&ione,w2,&ione);
//                 dcjacapply(jacs, n, rx1, rx2, w2, w);
//                 dbfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,w,w);
            }
            else {
                dbfun3(&rx1,&n,&rx2, &rx1,&n,&rx2, &ra1, &ra2, Phi1,A,Phi2, w, w, res1, res2);
//                 dbfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,w,w);
            }

            // Orthog w to j projectors
            for (i=0; i<=j; i++) {
                dbeta = -2.0*ddot_(&sz,&U[i*sz],&ione,w,&ione);
                daxpy_(&sz,&dbeta,&U[i*sz],&ione,w,&ione);
            }

            // new P_{j+1}
            if (j<sz-1) {
                for (i=0; i<=j; i++) U[i+(j+1)*sz]=0.0;
                i = sz-j-1;
                dcopy_(&i, &w[j+1], &ione, &U[j+1+(j+1)*sz], &ione);
                dalpha = dnrm2_(&i, &U[j+1+(j+1)*sz], &ione);
                if (dalpha!=0.0) {
                    if (w[j+1]<0.0) dalpha = -dalpha;
                    U[j+1+(j+1)*sz]+=dalpha;
                    dbeta = dnrm2_(&i, &U[j+1+(j+1)*sz], &ione);
                    dbeta = 1.0/dbeta;
                    dscal_(&i,&dbeta,&U[j+1+(j+1)*sz],&ione);

                    w[j+1]=-dalpha;
                    for (i=j+2; i<sz; i++) w[i]=0.0;
                }
            }

            // Givens rotators to the top of w
            for (i=0; i<=j-1; i++) {
                dbeta = w[i];
                w[i] = J[0+i*2]*w[i] + J[1+i*2]*w[i+1];
                w[i+1] = -J[1+i*2]*dbeta + J[0+i*2]*w[i+1];
            }

//            mexPrintf("w[%d]=%3.7e\tw[%d]=%3.7e\n", i, w[j], j+1, w[j+1]);

            // New rotator
            if (j<sz-1) {
                dalpha = sqrt((w[j]*w[j])+(w[j+1]*w[j+1]));
//                 mexPrintf("rho=%3.7e\n", dalpha);
                J[0+j*2] = w[j]/dalpha;
                J[1+j*2] = w[j+1]/dalpha;
//                 mexPrintf("Jnew=%3.7e\t%3.7e\n", J[j*2], J[j*2+1]);
                tau[j+1] = -J[1+j*2]*tau[j];
                tau[j] = J[0+j*2]*tau[j];
//                 mexPrintf("tau=%3.7e\t%3.7e\n", tau[j], tau[j+1]);
                w[j] = dalpha;
                w[j+1] = 0.0;
            }

//             for (i=0; i<=j+1; i++) mexPrintf("tau[%d]=%3.7e\n", i, tau[i]);

            dcopy_(&nrestart, w, &ione, &R[j*nrestart], &ione);


            // residual
            curres = fabs(tau[j+1])/nrmrhs;
            if (verb>1) printf("iter [%d,%d], res: %3.5e\n", it, j, curres);
//             if (verb>1) sprintf(&gmres_3d_printf_buf_[0], "iter [%d,%d], res: %3.5e\n", it, j, curres);

            if ((curres<tol)) break;
        }

//         for (i=0; i<=j; i++) mexPrintf("%g\t", tau[i]);
//         mexPrintf("\n");

        if (j==nrestart) {
            j=nrestart-1;
            i=nrestart;
        }
        else i=j+1;
        dtrsv_(&uplo,&trans,&trans,&i,R,&nrestart,tau,&ione);

/*        if (j==nrestart) {
            j=nrestart-1;
            i=nrestart;
        }
        else i=j+1;*/

        // Correction
        dcopy_(&sz, &U[j*sz], &ione, w, &ione);
        dbeta = -2.0*U[j+j*sz]*tau[j];
        dscal_(&sz, &dbeta, w, &ione);
        w[j]+=tau[j];
        for (i=j-1; i>=0; i--) {
            w[i]+=tau[i];
            dbeta = -2.0*ddot_(&sz,&U[i*sz],&ione,w,&ione);
            daxpy_(&sz,&dbeta,&U[i*sz],&ione,w,&ione);
        }
        dalpha=1.0;
        daxpy_(&sz,&dalpha,w,&ione,sol,&ione);
        if ((curres<tol)) break;
    }
//     if (verb>0) mexPrintf("gmres conducted [%d,%d] iters to relres %3.3e\n", it, j, curres);
    if (verb>0) sprintf(&gmres_3d_printf_buf_[0], "gmres conducted [%d,%d] iters to relres %3.3e\n", it, j, curres);
//     if (verb>0) printf("gmres conducted [%d,%d] iters to relres %3.3e\n", it, j, curres);

    free(U);
    free(w);
//     free(w2);
    free(tau);
    free(J);
    free(R);
    free(res1);
    free(res2);
}


void dgmresl_hh(double *Phi1, double *A, double *Phi2, double *rhs, long rx1, long n, long rx2, long ra1, long ra2, long nrestart, double tol, long niters, char ptype, double *jacs, double *sol, char verb)
// left preconditioned - for fro tolerance
// This one with Householder tranforms
{
    long i,j,it, sz;
    double *U, *w;
    double *R, *J, *tau;
    double *res1, *res2;
    double nrmr, curres, nrmrhs, nrmsol;
    double dalpha, dbeta;
    char last_iter = 0;
    long ione = 1;
    char trans = 'N';
    char uplo = 'U';

    sz = rx1*n*ra2*rx2;
    if (rx1*ra1*n*rx2>sz) sz = rx1*ra1*n*rx2;
    res1 = (double *)malloc(sizeof(double)*sz);
    res2 = (double *)malloc(sizeof(double)*sz);

    sz = rx1*n*rx2;

    U = (double *)malloc(sizeof(double)*sz*(nrestart+1));
    w = (double *)malloc(sizeof(double)*sz);
//     w2 = (double *)malloc(sizeof(double)*sz);
    tau = (double *)malloc(sizeof(double)*(nrestart+1));
    J = (double *)malloc(sizeof(double)*2*(nrestart));
//     Q = (double *)malloc(sizeof(double)*nrestart*(nrestart+1));
    R = (double *)malloc(sizeof(double)*(nrestart)*(nrestart));

    if ((ptype!='n')&&(ptype!='N')) {
        dcopy_(&sz,rhs,&ione,w,&ione);
        djac_apply(&ptype, &rx1, &n, &rx2, jacs, w, w, res1);
//         dcjacapply(jacs, n, rx1, rx2, w, w2);
        nrmrhs = dnrm2_(&sz, w, &ione);
    } else
        nrmrhs = dnrm2_(&sz, rhs, &ione);

    for (it=0; it<niters; it++) {
        // r0
        nrmsol = dnrm2_(&sz, sol, &ione);
        dbfun3(&rx1,&n,&rx2, &rx1,&n,&rx2, &ra1, &ra2, Phi1,A,Phi2, sol, w, res1, res2);
//         dbfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,sol,w);
        dbeta = -1.0;
        daxpy_(&sz,&dbeta,rhs,&ione,w,&ione);
        dscal_(&sz,&dbeta,w,&ione);
        if ((ptype!='n')&&(ptype!='N')) {
//             dcopy_(&sz,w,&ione,w2,&ione);
//             dcjacapply(jacs, n, rx1, rx2, w2, w);
            djac_apply(&ptype, &rx1, &n, &rx2, jacs, w, w, res1);
        }
        nrmr = dnrm2_(&sz,w,&ione);
        if (verb>1) printf("restart %d, res: %3.5e\n", it, nrmr);
//         if (verb>1) sprintf(&gmres_3d_printf_buf_[0], "restart %d, res: %3.5e\n", it, nrmr);
        if (it==0) nrmrhs = nrmr;
        // initial HHT
        dbeta = nrmr;
        if (w[0]<0.0) dbeta = -dbeta;
        w[0]+=dbeta;
        tau[0]=-dbeta;
        nrmr = dnrm2_(&sz,w,&ione);
        dbeta = 1.0/nrmr;
        dscal_(&sz,&dbeta,w,&ione);
        dcopy_(&sz,w,&ione,&U[0],&ione);

        for (j=0; j<nrestart; j++) {
            // HHT on last U
            dcopy_(&sz,&U[sz*j],&ione,w,&ione);
            dbeta = -2.0*U[j+sz*j];
            dscal_(&sz,&dbeta,w,&ione);
            w[j]+=1.0;
            for (i=j-1; i>=0; i--) {
                dbeta = -2.0*ddot_(&sz,&U[i*sz],&ione,w,&ione);
                daxpy_(&sz,&dbeta,&U[i*sz],&ione,w,&ione);
            }
            dbeta = dnrm2_(&sz,w,&ione);
            dbeta = 1.0/dbeta;
            dscal_(&sz,&dbeta,w,&ione); // w=w/norm(w);

//             mexPrintf("norm_v0=%3.7e\n", 1.0/dbeta);

            // precvec, matvec
            dbfun3(&rx1,&n,&rx2, &rx1,&n,&rx2, &ra1, &ra2, Phi1,A,Phi2, w, w, res1, res2);
//             dbfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,w,w);
            if ((ptype!='n')&&(ptype!='N')) {
                djac_apply(&ptype, &rx1, &n, &rx2, jacs, w, w, res1);
/*                dcopy_(&sz,w,&ione,w2,&ione);
                dcjacapply(jacs, n, rx1, rx2, w2, w);*/
            }

            // Orthog w to j projectors
            for (i=0; i<=j; i++) {
                dbeta = -2.0*ddot_(&sz,&U[i*sz],&ione,w,&ione);
                daxpy_(&sz,&dbeta,&U[i*sz],&ione,w,&ione);
            }

            // new P_{j+1}
            if (j<sz-1) {
                for (i=0; i<=j; i++) U[i+(j+1)*sz]=0.0;
                i = sz-j-1;
                dcopy_(&i, &w[j+1], &ione, &U[j+1+(j+1)*sz], &ione);
                dalpha = dnrm2_(&i, &U[j+1+(j+1)*sz], &ione);
                if (dalpha!=0.0) {
                    if (w[j+1]<0.0) dalpha = -dalpha;
                    U[j+1+(j+1)*sz]+=dalpha;
                    dbeta = dnrm2_(&i, &U[j+1+(j+1)*sz], &ione);
                    dbeta = 1.0/dbeta;
                    dscal_(&i,&dbeta,&U[j+1+(j+1)*sz],&ione);

                    w[j+1]=-dalpha;
                    for (i=j+2; i<sz; i++) w[i]=0.0;
                }
            }

            // Givens rotators to the top of w
            for (i=0; i<=j-1; i++) {
                dbeta = w[i];
                w[i] = J[0+i*2]*w[i] + J[1+i*2]*w[i+1];
                w[i+1] = -J[1+i*2]*dbeta + J[0+i*2]*w[i+1];
            }

//            mexPrintf("w[%d]=%3.7e\tw[%d]=%3.7e\n", i, w[j], j+1, w[j+1]);

            // New rotator
            if (j<sz-1) {
                dalpha = sqrt((w[j]*w[j])+(w[j+1]*w[j+1]));
//                 mexPrintf("rho=%3.7e\n", dalpha);
                J[0+j*2] = w[j]/dalpha;
                J[1+j*2] = w[j+1]/dalpha;
//                 mexPrintf("Jnew=%3.7e\t%3.7e\n", J[j*2], J[j*2+1]);
                tau[j+1] = -J[1+j*2]*tau[j];
                tau[j] = J[0+j*2]*tau[j];
//                 mexPrintf("tau=%3.7e\t%3.7e\n", tau[j], tau[j+1]);
                w[j] = dalpha;
                w[j+1] = 0.0;
            }

//             for (i=0; i<=j+1; i++) mexPrintf("tau[%d]=%3.7e\n", i, tau[i]);

            dcopy_(&nrestart, w, &ione, &R[j*nrestart], &ione);


            // residual
            curres = fabs(tau[j+1])/nrmrhs;
            if (verb>1) printf("iter [%d,%d], res: %3.5e\n", it, j, curres);
//             if (verb>1) sprintf(&gmres_3d_printf_buf_[0], "iter [%d,%d], res: %3.5e\n", it, j, curres);

            if ((curres<tol)) break;
        }

//         for (i=0; i<=j; i++) mexPrintf("%g\t", tau[i]);
//         mexPrintf("\n");

        if (j==nrestart) {
            j=nrestart-1;
            i=nrestart;
        }
        else i=j+1;
        dtrsv_(&uplo,&trans,&trans,&i,R,&nrestart,tau,&ione);

/*        if (j==nrestart) {
            j=nrestart-1;
            i=nrestart;
        }
        else i=j+1;*/

        // Correction
        dcopy_(&sz, &U[j*sz], &ione, w, &ione);
        dbeta = -2.0*U[j+j*sz]*tau[j];
        dscal_(&sz, &dbeta, w, &ione);
        w[j]+=tau[j];
        for (i=j-1; i>=0; i--) {
            w[i]+=tau[i];
            dbeta = -2.0*ddot_(&sz,&U[i*sz],&ione,w,&ione);
            daxpy_(&sz,&dbeta,&U[i*sz],&ione,w,&ione);
        }
        dalpha=1.0;
//         dcopy_(&sz, w, &ione, sol, &ione);
        daxpy_(&sz,&dalpha,w,&ione,sol,&ione);
        if ((curres<tol)&&(fabs(tau[j])/nrmsol<tol)) break;
    }
//     if (verb>0) mexPrintf("gmres conducted [%d,%d] iters to relres %3.3e\n", it, j, curres);
    if (verb>0) sprintf(&gmres_3d_printf_buf_[0], "gmres conducted [%d,%d] iters to relres %3.3e\n", it, j, curres);

    free(U);
    free(w);
//     free(w2);
    free(tau);
    free(J);
    free(R);
    free(res1);
    free(res2);
}


// Fortran wrappers
void dgmresr_hh_(double *Phi1, double *A, double *Phi2, double *rhs, long *rx1, long *n, long *rx2, long *ra1, long *ra2, long *nrestart, double *tol, long *niters, char *ptype, double *jacs, double *sol, char *verb)
{
//     printf("sz: %dx%dx%d, ra:%dx%d, its: %dx%d, prec: %c, verb:%d, tol:%g\n",rx1[0],n[0],rx2[0],ra1[0],ra2[0],nrestart[0],niters[0], ptype[0], verb[0], tol[0]);

    dgmresr_hh(Phi1, A, Phi2, rhs, rx1[0], n[0], rx2[0], ra1[0], ra2[0], nrestart[0], tol[0], niters[0], ptype[0], jacs, sol, verb[0]);
}

void dgmresl_hh_(double *Phi1, double *A, double *Phi2, double *rhs, long *rx1, long *n, long *rx2, long *ra1, long *ra2, long *nrestart, double *tol, long *niters, char *ptype, double *jacs, double *sol, char *verb)
{
    dgmresl_hh(Phi1, A, Phi2, rhs, rx1[0], n[0], rx2[0], ra1[0], ra2[0], nrestart[0], tol[0], niters[0], ptype[0], jacs, sol, verb[0]);
}
