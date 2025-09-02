#ifndef MATHIEU_EIGS_H
#define MATHIEU_EIGS_H

#include "../config.h"
#include "../error.h"
#include "make_matrix.h"
#include "matrix_utils.h"


/*
 *
 * This is part of the Mathieu function suite -- a reimplementation
 * of the Mathieu functions for Scipy.  This file holds the functions
 * which return the Mathieu eigenvalues (characteristic values) a and
 * b as a function of parameter q.
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
  int mathieu_a(int m, double q, double *a) {
    int N = m+25;  // Sets size of recursion matrix
    int retcode = SF_ERROR_OK;

    // Allocate recursion matrix
    double *A = (double *) calloc(N*N, sizeof(double));
    if (A == NULL) return SF_ERROR_MEMORY;
 
    // Allocate vector for eigenvalues
    double *ww = (double *) calloc(N, sizeof(double));
    if (ww == NULL) {
      free(A);
      return SF_ERROR_MEMORY;
    }

    // Do EVD
    if (m % 2 == 0) {
      // Even order m
      retcode = make_matrix_ee(N,q,A);
      if (retcode != 0){
	free(A);
	free(ww);
	return retcode;
      }
      {
	char V = 'V';
	char U = 'U';      
	double wkopt;
	double* work;
	/* Query and allocate the optimal workspace */
	// I do this work in an inner scope to make it easy
	// to clean up afterwards.
	int lwork = -1;
	dsyev_( &V, &U, &N, A, &N, ww, &wkopt, &lwork, &retcode );
	lwork = (int) wkopt;
	work = (double*)malloc( lwork*sizeof(double) );
	/* Solve eigenproblem */
	dsyev_( &V, &U, &N, A, &N, ww, work, &lwork, &retcode );
	free(work);
      }
      if (retcode != 0) {
	free(A);
	free(ww);
	return SF_ERROR_NO_RESULT;
      }
      
      // Sort ww vector from lowest to highest
      quickSort(ww, 0, N-1);
      //print_matrix(ww, N, 1);      

      // Now figure out which one to return.
      int idx = m/2;
      *a = ww[idx];
      
    } else {
      // Odd order m
      retcode = make_matrix_eo(N,q,A);
      if (retcode != 0) {
	free(A);
	free(ww);
	return retcode;
      }
      {
	char V = 'V';
	char U = 'U';      
	double wkopt;
	double* work;
	/* Query and allocate the optimal workspace */
	// I do this work in an inner scope to make it easy
	// to clean up afterwards.
	int lwork = -1;
	dsyev_( &V, &U, &N, A, &N, ww, &wkopt, &lwork, &retcode );
	lwork = (int) wkopt;
	work = (double*)malloc( lwork*sizeof(double) );
	/* Solve eigenproblem */
	dsyev_( &V, &U, &N, A, &N, ww, work, &lwork, &retcode );
	free(work);
      }
      if (retcode != 0) {
	free(A);
	free(ww);
	return SF_ERROR_NO_RESULT;
      }
      
      // Sort ww vector from lowest to highest
      quickSort(ww, 0, N-1);
      //print_matrix(ww, N, 1);
      
      // Now figure out which one to return.
      int idx = (m-1)/2;
      *a = ww[idx];
    }

    free(A);
    free(ww);
    
    return retcode;
  }

  //------------------------------------------------------
  int mathieu_b(int m, double q, double *b) {
    int N = m+25;  // Sets size of recursion matrix
    int retcode = SF_ERROR_OK;

    // Allocate recursion matrix
    double *B = (double *) calloc(N*N, sizeof(double));
    if (B == NULL) return SF_ERROR_MEMORY;

    // Allocate vector for eigenvalues
    double *ww = (double *) calloc(N, sizeof(double));
    if (ww == NULL) {
      free(B);
      return SF_ERROR_MEMORY;
    }

    // Do EVD
    if (m % 2 == 0) {
      // Even order m
      retcode = make_matrix_oe(N,q,B);
      if (retcode != 0) {
	free(B);
	free(ww);
	return retcode;
      }
      {
	char V = 'V';
	char U = 'U';      
	double wkopt;
	double* work;
 	/* Query and allocate the optimal workspace */
	// I do this work in an inner scope to make it easy
	// to clean up afterwards.
	int lwork = -1;
	dsyev_( &V, &U, &N, B, &N, ww, &wkopt, &lwork, &retcode );
	lwork = (int) wkopt;
	work = (double*)malloc( lwork*sizeof(double) );
	/* Solve eigenproblem */
	dsyev_( &V, &U, &N, B, &N, ww, work, &lwork, &retcode );
	free(work);
      }
      if (retcode != 0) {
	free(B);
	free(ww);
	return SF_ERROR_NO_RESULT;
      }
      
      // Sort ww vector from lowest to highest
      quickSort(ww, 0, N-1);
      //print_matrix(ww, N, 1);
      
      // Now figure out which one to return.
      int idx = (m-2)/2;
      *b = ww[idx];
      
    } else {
      // Odd order m
      retcode = make_matrix_oo(N,q,B);
      if (retcode != 0) {
	free(B);
	free(ww);
	return retcode;
      }
      {
	char V = 'V';
	char U = 'U';      
	double wkopt;
	double* work;
	/* Query and allocate the optimal workspace */
	// I do this work in an inner scope to make it easy
	// to clean up afterwards.
	int lwork = -1;
	dsyev_( &V, &U, &N, B, &N, ww, &wkopt, &lwork, &retcode );
	lwork = (int) wkopt;
	work = (double*)malloc( lwork*sizeof(double) );
	/* Solve eigenproblem */
	dsyev_( &V, &U, &N, B, &N, ww, work, &lwork, &retcode );
	free(work);
      }
      if (retcode != 0) {
	free(B);
	free(ww);
	return SF_ERROR_NO_RESULT;
      }
      
      // Sort ww vector from lowest to highest
      quickSort(ww, 0, N-1);
      //print_matrix(ww, N, 1);

      // Now figure out which one to return.
      int idx = (m-1)/2;
      *b = ww[idx];
      
    }

    free(B);
    free(ww);
    return retcode;
  }
    
} // namespace mathieu
} // namespace xsf

#endif // #ifndef MATHIEU_EIGS_H
