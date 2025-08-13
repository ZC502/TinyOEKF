/*
 * Octonion Extended Kalman Filter for embedded processors
 *
 * Copyright (C) 2024 Simon D. Levy
 * Modifications Copyright (C) 2025 ZuoCen Liu
 * MIT License
 *
 */

#include <math.h>
#include <stdbool.h>
#include <string.h>

/**
  * Floating-point precision defaults to single but can be made double via
    <tt><b>#define _float_t double</b></tt> before <tt>#include <tinyekf.h></tt>
  */
#ifndef _float_t
#define _float_t float
#endif

// define OEKF dimension macros
#define OEKF_N 14     // Total state dimensions: 8 (octonions) + 3 (velocity) + 3 (position)
#define OEKF_M 6    // Observation dimensions (defined according to the actual number of sensors, e.g., 6 for IMU + GPS)

// Linear alegbra ////////////////////////////////////////////////////////////

/// @private
static void _mulmat(
        const _float_t * a, 
        const _float_t * b, 
        _float_t * c, 
        const int arows, 
        const int acols, 
        const int bcols)
{
    for (int i=0; i<arows; ++i) {
        for (int j=0; j<bcols; ++j) {
            c[i*bcols+j] = 0;
            for (int k=0; k<acols; ++k) {
                c[i*bcols+j] += a[i*acols+k] * b[k*bcols+j];
            }
        }
    }
}

/// @private
static void _mulvec(
        const _float_t * a, 
        const _float_t * x, 
        _float_t * y, 
        const int m, 
        const int n)
{
    for (int i=0; i<m; ++i) {
        y[i] = 0;
        for (int j=0; j<n; ++j)
            y[i] += x[j] * a[i*n+j];
    }
}

/// @private
static void _transpose(
        const _float_t * a, _float_t * at, const int m, const int n)
{
    for (int i=0; i<m; ++i)
        for (int j=0; j<n; ++j) {
            at[j*m+i] = a[i*n+j];
        }
}

/// @private
static void _addmat(
        const _float_t * a, const _float_t * b, _float_t * c, 
        const int m, const int n)
{
    for (int i=0; i<m; ++i) {
        for (int j=0; j<n; ++j) {
            c[i*n+j] = a[i*n+j] + b[i*n+j];
        }
    }
}

/// @private
static void _negate(_float_t * a, const int m, const int n)
{        
    for (int i=0; i<m; ++i) {
        for (int j=0; j<n; ++j) {
            a[i*n+j] = -a[i*n+j];
        }
    }
}

/// @private
static void _addeye(_float_t * a, const int n)
{
    for (int i=0; i<n; ++i) {
        a[i*n+i] += 1;
    }
}


/* Cholesky-decomposition matrix-inversion code, adapated from
http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/_choles_cpp.txt */

/// @private
static int _choldc1(_float_t * a, _float_t * p, const int n) 
{
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            _float_t sum = a[i*n+j];
            for (int k = i - 1; k >= 0; k--) {
                sum -= a[i*n+k] * a[j*n+k];
            }
            if (i == j) {
                if (sum <= 0) {
                    return 1; /* error */
                }
                p[i] = sqrt(sum);
            }
            else {
                a[j*n+i] = sum / p[i];
            }
        }
    }

    return 0; // success:w
}

/// @private
static int _choldcsl(const _float_t * A, _float_t * a, _float_t * p, const int n) 
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i*n+j] = A[i*n+j];
        }
    }
    if (_choldc1(a, p, n)) {
        return 1;
    }
    for (int i = 0; i < n; i++) {
        a[i*n+i] = 1 / p[i];
        for (int j = i + 1; j < n; j++) {
            _float_t sum = 0;
            for (int k = i; k < j; k++) {
                sum -= a[j*n+k] * a[k*n+i];
            }
            a[j*n+i] = sum / p[j];
        }
    }

    return 0; // success
}

/// @private
static int _cholsl(const _float_t * A, _float_t * a, _float_t * p, const int n) 
{
    if (_choldcsl(A,a,p,n)) {
        return 1;
    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            a[i*n+j] = 0.0;
        }
    }
    for (int i = 0; i < n; i++) {
        a[i*n+i] *= a[i*n+i];
        for (int k = i + 1; k < n; k++) {
            a[i*n+i] += a[k*n+i] * a[k*n+i];
        }
        for (int j = i + 1; j < n; j++) {
            for (int k = j; k < n; k++) {
                a[i*n+j] += a[k*n+i] * a[k*n+j];
            }
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            a[i*n+j] = a[j*n+i];
        }
    }

    return 0; // success
}

/// @private
static void _addvec(
        const _float_t * a, const _float_t * b, _float_t * c, const int n)
{
    for (int j=0; j<n; ++j) {
        c[j] = a[j] + b[j];
    }
}

/// @private
static void _sub(
        const _float_t * a, const _float_t * b, _float_t * c, const int n)
{
    for (int j=0; j<n; ++j) {
        c[j] = a[j] - b[j];
    }
}

/// @private
static bool invert(const _float_t * a, _float_t * ainv)
{
    _float_t tmp[OEKF_M];

    return _cholsl(a, ainv, tmp, OEKF_M) == 0;
}


// OEKF ///////////////////////////////////////////////////////////////////////

// New: Octonion structure (1 real part + 7 imaginary parts, single precision)
typedef struct {
 _float_t r;        // Real part (same precision as the original code)
_float_t i[7];      // Imaginary parts (i0~i6, corresponding to 7 dimensions)
} Octonion;

// New: OEKF state structure (octonion + velocity + position, total dimensions = 8 + 3 + 3 = 14)
typedef struct {
Octonion q;        // Octonion (rotation + translation fusion, 8D)
_float_t v[3];     // Velocity (x, y, z, 3D)
_float_t p[3];     // Position (x, y, z, 3D)
} OEKF_State;

// Replace the original ekf_t with a structure dedicated to OKF
typedef struct {
OEKF_State state;        // State (quaternion + velocity + position)
_float_t P[14 * 14];    // Covariance matrix (14×14, matching the state dimension)
}oekf_t;

static void okf_initialize(okf_t *okf, const _float_t pdiag[OKF_N]) {  
 //  Initialize the octonion to the identity element (no rotation)
    okf->state.q.r = 1.0f;
    memset(okf->state.q.i, 0, 7*sizeof(_float_t));
    // Initialize the velocity and position to 0
    memset(okf->state.v, 0, 3*sizeof(_float_t));
    memset(okf->state.p, 0, 3*sizeof(_float_t));
    // Initialize the covariance matrix (diagonal matrix)
    for (int i=0; i<OKF_N; ++i) {
        for (int j=0; j<OKF_N; ++j) {
            okf->P[i*OKF_N + j] = (i==j) ? pdiag[i] : 0;
        }
    }
}

/**
 * Initializes the EKF
 * @param ekf pointer to an oekf_t structure
 * @param pdiag a vector of length OEKF_N containing the initial values for the
 * covariance matrix diagonal
 */


/**
  * Runs the EKF prediction step
  * @param ekf pointer to an oekf_t structure
  * @param fx predicted values
  * @param F Jacobian of state-transition function
  * @param Q process noise matrix
  * 
  */static void oekf_predict(
       oekf_t * ekf, 
        const _float_t fx[OEKF_N],
        const _float_t F[OEKF_N*OEKF_N],
        const _float_t Q[OEKF_N*OEKF_N])
{        
    // \hat{x}_k = f(\hat{x}_{k-1}, u_k)
    memcpy(ekf->x, fx, OEKF_N*sizeof(_float_t));

    // P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1}

    _float_t FP[OEKF_N*OEKF_N] = {};
    _mulmat(F, ekf->P,  FP, OEKF_N, OEKF_N, OEKF_N);

    _float_t Ft[OEKF_N*OEKF_N] = {};
    _transpose(F, Ft, OEKF_N, OEKF_N);

    _float_t FPFt[OEKF_N*OEKF_N] = {};
    _mulmat(FP, Ft, FPFt, OEKF_N, OEKF_N, OEKF_N);

    _addmat(FPFt, Q, ekf->P, OEKF_N, OEKF_N);
}

/// @private
static void oekf_update_step3(oekf_t * ekf, _float_t GH[OEKF_N*OEKF_N])
{
    _negate(GH, OEKF_N, OEKF_N);
    _addeye(GH, OEKF_N);
    _float_t GHP[OEKF_N*OEKF_N];
    _mulmat(GH, ekf->P, GHP, OEKF_N, OEKF_N, OEKF_N);
    memcpy(ekf->P, GHP, OEKF_N*OEKF_N*sizeof(_float_t));
}

/**
  * Runs the EKF update step
  * @param ekf pointer to an oekf_t structure
  * @param z observations
  * @param hx predicted values
  * @param H sensor-function Jacobian matrix
  * @param R measurement-noise matrix
  * 
  */
static bool oekf_update(
       oekf_t * ekf, 
        const _float_t z[OEKF_M], 
        const _float_t hx[OEKF_N],
        const _float_t H[OEKF_M*OEKF_N],
        const _float_t R[OEKF_M*OEKF_M])
{        
    // G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}
    _float_t G[OEKF_N*OEKF_M];
    _float_t Ht[OEKF_N*OEKF_M];
    _transpose(H, Ht, OEKF_M, OEKF_N);
    _float_t PHt[OEKF_N*OEKF_M];
    _mulmat(ekf->P, Ht, PHt, OEKF_N, OEKF_N, OEKF_M);
    _float_t HP[OEKF_M*OEKF_N];
    _mulmat(H, ekf->P, HP, OEKF_M, OEKF_N, OEKF_N);
    _float_t HpHt[OEKF_M*OEKF_M];
    _mulmat(HP, Ht, HpHt, OEKF_M, OEKF_N, OEKF_M);
    _float_t HpHtR[OEKF_M*OEKF_M];
    _addmat(HpHt, R, HpHtR, OEKF_M, OEKF_M);
    _float_t HPHtRinv[OEKF_M*OEKF_M];
    if (!invert(HpHtR, HPHtRinv)) {
        return false;
    }
    _mulmat(PHt, HPHtRinv, G, OEKF_N, OEKF_M, OEKF_M);

    // \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))
    _float_t z_hx[OEKF_M];
    _sub(z, hx, z_hx, OEKF_M);
    _float_t Gz_hx[OEKF_M*OEKF_N];
    _mulvec(G, z_hx, Gz_hx, OEKF_N, OEKF_M);
    _addvec(ekf->x, Gz_hx, ekf->x, OEKF_N);

    // P_k = (I - G_k H_k) P_k
    _float_t GH[OEKF_N*OEKF_N];
    _mulmat(G, H, GH, OEKF_N, OEKF_M, OEKF_N);
    oekf_update_step3(ekf, GH);

    // success
    return true;
}
