#ifndef MATHIEU_COEFFS_H
#define MATHIEU_COEFFS_H

#include "../config.h"
#include "../error.h"
#include "make_matrix.h"
#include "matrix_utils.h"
#include "mathieu_coeffs.h"
// #include <cblas.h>
// #include <lapacke.h>
// #include <lapack.h>

#define SQRT2 1.414213562373095d

/*
 *
 * This is part of the Mathieu function suite -- a reimplementation
 * of the Mathieu functions for Scipy.  This file holds the functions
 * which return the Mathieu A and B coefficients used in the Fourier
 * series computing the Mathieu fcns.
 * 
 */



/* DSYEV_ prototype */
#ifdef __cplusplus
extern "C" {
#endif
void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,
	    double* w, double* work, int* lwork, int* info );
#ifdef __cplusplus
}
#endif

  
namespace xsf {
namespace mathieu {

  //------------------------------------------------------
  int mathieu_coeffs_ee(int N, double q, int m, double *AA) {
    // Returns Fourier coeffs for the mth order ce_2n Mathieu fcn.
    // Allowed value of m = 0, 2, 4, 6, ...
    // Inputs:
    // N = size of recursion matrix to use.
    // q = frequency parameter
    // m = order of Mathieu fcn desired.
    // Output:
    // AA = length N vector preallocated to hold coeffs.
    // Returns 0 if all goes well.  Must put check on calloc
    // here.

    int retcode = 0;

    // Bail out if m is not even.
    if (m % 2 != 0) return -1;

    // Allocate recursion matrix
    double *A = (double *) calloc(N*N, sizeof(double));
    if (A == NULL) return SF_ERROR_MEMORY;

    // Do EVD
    //printf("Even m = %d\n", m);
    retcode = make_matrix_ee(N,q,A);
    if (retcode != 0) {
      free(A);
      return SF_ERROR_NO_RESULT;
    }
    //print_matrix(A, N, N);
    //printf("---------------------------------------\n");
    {
      char V[1] = {'V'};
      char U[1] = {'U'};      
      double wkopt;
      double* work;
      /* Query and allocate the optimal workspace */
      // I do this work in an inner scope to make it easy
      // to clean up afterwards.
      int lwork = -1;
      dsyev_( V, U, &N, A, &N, AA, &wkopt, &lwork, &retcode );
      lwork = (int)wkopt;
      work = (double*)malloc( lwork*sizeof(double) );
      /* Solve eigenproblem */
      dsyev_( V, U, &N, A, &N, AA, work, &lwork, &retcode );
      free(work);
    }
    // I replaced LAPACKE call with above LAPACK call to integrate into Scipy.
    // retcode = LAPACKE_dsyev_(LAPACK_ROW_MAJOR, 'V', 'L', N, A, N, AA);
    if (retcode != 0) {
	free(A);
	return SF_ERROR_NO_RESULT;
    }
    //print_matrix(A, N, N);
    //printf("---------------------------------------\n");    
    
    // Sort AA vector from lowest to highest
    // quickSort(AA, 0, N-1);
    //print_matrix(AA, N, 1);
    
    // Undo sqrt(2) in make_matrix by normalizing elet in first col by sqrt(2).
    int idx;
    int row = m/2;
    idx = MATRIX_IDX(N, row, 0);
    AA[0] = A[idx]/SQRT2;
    // Transfer remaining elets in correct row to coeff vector.
    for (int j = 1; j < N; j++) {
      idx = MATRIX_IDX(N, row, j);
      AA[j] = A[idx];
    }
    free(A);
    return retcode;
  }


  //------------------------------------------------------
  int mathieu_coeffs_eo(int N, double q, int m, double *AA) {
    // Returns Fourier coeffs for the mth order ce_2n+1 Mathieu fcn.
    // Allowed value of m = 1, 3, 5, 7 ... 

    int retcode = 0;

    // Bail out if m is not odd.
    if (m % 2 != 1) return -1;

    // Allocate recursion matrix
    double *A = (double *) calloc(N*N, sizeof(double));
    if (A == NULL) return SF_ERROR_MEMORY;
    
    // Do EVD
    //printf("Odd m = %d\n", m);
    retcode = make_matrix_eo(N,q,A);
    if (retcode != 0) {
      free(A);
      return SF_ERROR_NO_RESULT;
    }
    //print_matrix(A, N, N);
    //printf("---------------------------------------\n");
    {
      char V[1] = {'V'};
      char U[1] = {'U'};      
      double wkopt;
      double* work;
      /* Query and allocate the optimal workspace */
      // I do this work in an inner scope to make it easy
      // to clean up afterwards.
      int lwork = -1;
      dsyev_( V, U, &N, A, &N, AA, &wkopt, &lwork, &retcode );
      lwork = (int)wkopt;
      work = (double*)malloc( lwork*sizeof(double) );
      /* Solve eigenproblem */
      dsyev_( V, U, &N, A, &N, AA, work, &lwork, &retcode );
      free(work);
    }
    // I replaced LAPACKE call with above LAPACK call to integrate into Scipy.
    // retcode = LAPACKE_dsyev_(LAPACK_ROW_MAJOR, 'V', 'U', N, A, N, AA);
    if (retcode != 0) {
      free(A);
      return SF_ERROR_NO_RESULT;
    }
    //print_matrix(A, N, N);
    //printf("---------------------------------------\n");    
    
    // Sort AA vector from lowest to highest
    // quickSort(AA, 0, N-1);
    //print_matrix(AA, N, 1);

    // Transfer correct row to coeff vector.
    int idx;
    int row = (m-1)/2;
    // Transfer elets in correct row to coeff vector.
    for (int j = 0; j < N; j++) {
      idx = MATRIX_IDX(N, row, j);
      AA[j] = A[idx];
    }
    free(A);
    return retcode;
  }


  //------------------------------------------------------
  int mathieu_coeffs_oe(int N, double q, int m, double *AA) {
    // Returns Fourier coeffs for the mth order se_2n Mathieu fcn.
    // Allowed value of m = 2, 4, 6, ...
    // Inputs:
    // N = size of recursion matrix to use.
    // q = frequency parameter
    // m = order of Mathieu fcn desired.
    // Output:
    // AA = length N vector preallocated to hold coeffs.
    // Returns 0 if all goes well.  Must put check on calloc
    // here.

    int retcode = 0;

    // Bail out if m is not even or >= 2.
    if ((m % 2 != 0) || (m < 2)) return -1;

    // Allocate recursion matrix
    double *A = (double *) calloc(N*N, sizeof(double));
    if (A == NULL) return SF_ERROR_MEMORY;

    // Do EVD
    //printf("Even m = %d\n", m);
    retcode = make_matrix_oe(N,q,A);
    if (retcode != 0) {
      free(A);
      return SF_ERROR_NO_RESULT;
    }
    //print_matrix(A, N, N);
    //printf("---------------------------------------\n");
    {
      char V[1] = {'V'};
      char U[1] = {'U'};      
      double wkopt;
      double* work;
      /* Query and allocate the optimal workspace */
      // I do this work in an inner scope to make it easy
      // to clean up afterwards.
      int lwork = -1;
      dsyev_( V, U, &N, A, &N, AA, &wkopt, &lwork, &retcode );
      lwork = (int)wkopt;
      work = (double*)malloc( lwork*sizeof(double) );
      /* Solve eigenproblem */
      dsyev_( V, U, &N, A, &N, AA, work, &lwork, &retcode );
      free(work);
    }
    // I replaced LAPACKE call with above LAPACK call to integrate into Scipy.
    // retcode = LAPACKE_dsyev_(LAPACK_ROW_MAJOR, 'V', 'L', N, A, N, AA);
    if (retcode != 0) {
      free(A);
      return SF_ERROR_NO_RESULT;
    }
    //print_matrix(A, N, N);
    //printf("---------------------------------------\n");    
    
    // Sort AA vector from lowest to highest
    // quickSort(AA, 0, N-1);
    //print_matrix(AA, N, 1);
    
    // Transfer remaining elets in correct row to coeff vector.
    int idx;
    int row = (m-2)/2;
    for (int j = 0; j < N; j++) {
      idx = MATRIX_IDX(N, row, j);
      AA[j] = A[idx];
    }
    free(A);
    return retcode;
  }


  //------------------------------------------------------
  int mathieu_coeffs_oo(int N, double q, int m, double *AA) {
    // Returns Fourier coeffs for the mth order se_2n+1 Mathieu fcn.
    // Allowed value of m = 1, 3, 5, 7 ... 

    int retcode = 0;

    // Bail out if m is not odd.
    if (m % 2 != 1) return -1;

    // Allocate recursion matrix
    double *A = (double *) calloc(N*N, sizeof(double));
    if (A == NULL) return SF_ERROR_MEMORY;
    
    // Do EVD
    //printf("Odd m = %d\n", m);
    retcode = make_matrix_oo(N,q,A);
    if (retcode != 0) {
      free(A);
      return SF_ERROR_NO_RESULT;
    }
    //print_matrix(A, N, N);
    //printf("---------------------------------------\n");
    {
      char V[1] = {'V'};
      char U[1] = {'U'};      
      double wkopt;
      double* work;
      /* Query and allocate the optimal workspace */
      // I do this work in an inner scope to make it easy
      // to clean up afterwards.
      int lwork = -1;
      dsyev_( V, U, &N, A, &N, AA, &wkopt, &lwork, &retcode );
      lwork = (int)wkopt;
      work = (double*)malloc( lwork*sizeof(double) );
      /* Solve eigenproblem */
      dsyev_( V, U, &N, A, &N, AA, work, &lwork, &retcode );
      free(work);
    }
    // I replaced LAPACKE call with above LAPACK call to integrate into Scipy.
    // retcode = LAPACKE_dsyev_(LAPACK_ROW_MAJOR, 'V', 'U', N, A, N, AA);
    if (retcode != 0) {
      free(A);
      return SF_ERROR_NO_RESULT;
    }
    //print_matrix(A, N, N);
    //printf("---------------------------------------\n");    
    
    // Sort AA vector from lowest to highest
    // quickSort(AA, 0, N-1);
    //print_matrix(AA, N, 1);

    // Transfer correct row to coeff vector.
    int idx;
    int row = (m-1)/2;
    // Transfer elets in correct row to coeff vector.
    for (int j = 0; j < N; j++) {
      idx = MATRIX_IDX(N, row, j);
      AA[j] = A[idx];
    }
    free(A);
    return retcode;
  }


} // namespace mathieu
} // namespace xsf

//#ifdef __cplusplus
//}
//#endif
 
#endif
