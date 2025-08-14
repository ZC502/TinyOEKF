/*
 * Octonion Extended Kalman Filter for embedded processors
 *
 * Copyright (C) 2024 Simon D. Levy
 * Modifications Copyright (C) 2025 ZuoCen Liu
 * MIT License
 *
 */
#ifndef TINYOEKF_H
#define TINYOEKF_H

#include <math.h>
#include <stdbool.h>
#include <string.h>

#ifndef _float_t
#define _float_t float
#endif

#define OEKF_N 14  // 8(octonions) + 3(velocity) + 3(position)
#define OEKF_M 6   // Observation dimension

// Linear algebra tool functions
static void _mulmat(const _float_t *a, const _float_t *b, _float_t *c, const int arows, const int acols, const int bcols);
static void _mulvec(const _float_t *a, const _float_t *x, _float_t *y, const int m, const int n);
static void _transpose(const _float_t *a, _float_t *at, const int m, const int n);
static void _addmat(const _float_t *a, const _float_t *b, _float_t *c, const int m, const int n);
static void _negate(_float_t *a, const int m, const int n);
static void _addeye(_float_t *a, const int n);
static int _choldc1(_float_t *a, _float_t *p, const int n);
static int _choldcsl(const _float_t *A, _float_t *a, _float_t *p, const int n);
static int _cholsl(const _float_t *A, _float_t *a, _float_t *p, const int n);
static void _addvec(const _float_t *a, const _float_t *b, _float_t *c, const int n);
static void _sub(const _float_t *a, const _float_t *b, _float_t *c, const int n);
static bool invert(const _float_t *a, _float_t *ainv);

// Octonion structure
typedef struct {
    _float_t r;        // real part
    _float_t i[7];     // Imaginary part (i0~i6)
} Octonion;

// OEKF state structure
typedef struct {
    Octonion q;        // Octonions (8-dimensional)
    _float_t v[3];     // Velocity (3D)
    _float_t p[3];     // Position (3D)
} OEKF_State;

// Main structure of OEKF
typedef struct {
    OEKF_State state;  // State structure (for easy intuitive access)
    _float_t x[OEKF_N]; // State vector (used for EKF matrix operations)
    _float_t P[OEKF_N * OEKF_N]; // covariance matrix
} oekf_t;

// Initialization function
static void oekf_initialize(oekf_t *oekf, const _float_t pdiag[OEKF_N]) {
    // Initialize the octonion as the identity element
    oekf->state.q.r = 1.0f;
    memset(oekf->state.q.i, 0, 7 * sizeof(_float_t));
    
    // Initialize the velocity and position to 0
    memset(oekf->state.v, 0, 3 * sizeof(_float_t));
    memset(oekf->state.p, 0, 3 * sizeof(_float_t));
    
    // Map the state structure to the state vector x
    oekf->x[0] = oekf->state.q.r;
    memcpy(&oekf->x[1], oekf->state.q.i, 7 * sizeof(_float_t));
    memcpy(&oekf->x[8], oekf->state.v, 3 * sizeof(_float_t));
    memcpy(&oekf->x[11], oekf->state.p, 3 * sizeof(_float_t));
    
    // Initialize the covariance matrix (diagonal matrix)
    for (int i = 0; i < OEKF_N; ++i) {
        for (int j = 0; j < OEKF_N; ++j) {
            oekf->P[i * OEKF_N + j] = (i == j) ? pdiag[i] : 0;
        }
    }
}

// Prediction steps
static void oekf_predict(oekf_t *ekf, const _float_t fx[OEKF_N], const _float_t F[OEKF_N * OEKF_N], const _float_t Q[OEKF_N * OEKF_N]) {
    // Update the state vector
    memcpy(ekf->x, fx, OEKF_N * sizeof(_float_t));
    
    // Update the state structure (synchronized with vector x)
    ekf->state.q.r = ekf->x[0];
    memcpy(ekf->state.q.i, &ekf->x[1], 7 * sizeof(_float_t));
    memcpy(ekf->state.v, &ekf->x[8], 3 * sizeof(_float_t));
    memcpy(ekf->state.p, &ekf->x[11], 3 * sizeof(_float_t));
    
    // Calculate the covariance matrix P = F * P * F^T + Q
    _float_t FP[OEKF_N * OEKF_N] = {};
    _mulmat(F, ekf->P, FP, OEKF_N, OEKF_N, OEKF_N);
    
    _float_t Ft[OEKF_N * OEKF_N] = {};
    _transpose(F, Ft, OEKF_N, OEKF_N);
    
    _float_t FPFt[OEKF_N * OEKF_N] = {};
    _mulmat(FP, Ft, FPFt, OEKF_N, OEKF_N, OEKF_N);
    
    _addmat(FPFt, Q, ekf->P, OEKF_N, OEKF_N);
}

// Update step helper function
static void oekf_update_step3(oekf_t *ekf, _float_t GH[OEKF_N * OEKF_N]) {
    _negate(GH, OEKF_N, OEKF_N);
    _addeye(GH, OEKF_N);
    _float_t GHP[OEKF_N * OEKF_N];
    _mulmat(GH, ekf->P, GHP, OEKF_N, OEKF_N, OEKF_N);
    memcpy(ekf->P, GHP, OEKF_N * OEKF_N * sizeof(_float_t));
}

// Update steps
static bool oekf_update(oekf_t *ekf, const _float_t z[OEKF_M], const _float_t hx[OEKF_M], const _float_t H[OEKF_M * OEKF_N], const _float_t R[OEKF_M * OEKF_M]) {
    // Calculate the gain G = P * H^T * (H * P * H^T + R)^-1
    _float_t G[OEKF_N * OEKF_M];
    _float_t Ht[OEKF_N * OEKF_M];
    _transpose(H, Ht, OEKF_M, OEKF_N);
    
    _float_t PHt[OEKF_N * OEKF_M];
    _mulmat(ekf->P, Ht, PHt, OEKF_N, OEKF_N, OEKF_M);
    
    _float_t HP[OEKF_M * OEKF_N];
    _mulmat(H, ekf->P, HP, OEKF_M, OEKF_N, OEKF_N);
    
    _float_t HpHt[OEKF_M * OEKF_M];
    _mulmat(HP, Ht, HpHt, OEKF_M, OEKF_N, OEKF_M);
    
    _float_t HpHtR[OEKF_M * OEKF_M];
    _addmat(HpHt, R, HpHtR, OEKF_M, OEKF_M);
    
    _float_t HPHtRinv[OEKF_M * OEKF_M];
    if (!invert(HpHtR, HPHtRinv)) {
        return false;
    }
    
    _mulmat(PHt, HPHtRinv, G, OEKF_N, OEKF_M, OEKF_M);
    
    // Update the state vector x = x + G*(z - hx)
    _float_t z_hx[OEKF_M];
    _sub(z, hx, z_hx, OEKF_M);
    
    _float_t Gz_hx[OEKF_N];
    _mulvec(G, z_hx, Gz_hx, OEKF_N, OEKF_M);
    
    _addvec(ekf->x, Gz_hx, ekf->x, OEKF_N);
    
    // Update the state structure (synchronized with vector x)
    ekf->state.q.r = ekf->x[0];
    memcpy(ekf->state.q.i, &ekf->x[1], 7 * sizeof(_float_t));
    memcpy(ekf->state.v, &ekf->x[8], 3 * sizeof(_float_t));
    memcpy(ekf->state.p, &ekf->x[11], 3 * sizeof(_float_t));
    
    // Update the covariance matrix P = (I - G*H) * P
    _float_t GH[OEKF_N * OEKF_N];
    _mulmat(G, H, GH, OEKF_N, OEKF_M, OEKF_N);
    oekf_update_step3(ekf, GH);
    
    return true;
}

#endif // TINYOEKF_H

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
    _float_t tmp[EKF_M];

    return _cholsl(a, ainv, tmp, EKF_M) == 0;
}


// OEKF ///////////////////////////////////////////////////////////////////////

//Octonion structure (replacing the original quaternion/state vector)
typedef struct {
    float r;          // real part
    float i[7];       // 7 imaginary parts (i0-i6 correspond to e1-e7)
} Octonion;

// Octonion operation function
Octonion oct_multiply(const Octonion *a, const Octonion *b);
void oct_normalize(Octonion *q);

// Octonion Kalman filter state structure
typedef struct {
    Octonion q;       // Octonion state (rotation + translation)
    float P[14][14];  // 14-dimensional covariance matrix (8 quaternions + 3 velocities + 3 positions)
    float Q[14][14];  // Process noise matrix
    float R[6][6];    // Observation noise matrix (assuming 6-dimensional observation: 3 positions + 3 attitudes)
    float K[14][6];   // Kalman gain
    float dt;         // time step
} TinyOEKF;

typedef struct {
    OEKF_State state;                // Status structure
    _float_t x[OEKF_N];              // state vector
    _float_t P[OEKF_N * OEKF_N];     // covariance matrix
} oekf_t;

// Core function of the filter
void tinyoekf_init(TinyOEKF *ekf, float dt);
void tinyoekf_predict(TinyOEKF *ekf, const float *u);
void tinyoekf_update(TinyOEKF *ekf, const float *z);

/**
 * Initializes the EKF
 * @param ekf pointer to an ekf_t structure
 * @param pdiag a vector of length OEKF_N containing the initial values for the
 * covariance matrix diagonal
 */

// Synchronize state and x during initialization
static void oekf_initialize(oekf_t * oekf, const _float_t pdiag[OEKF_N]){
// Initialize the octonion as the identity element
   oekf->state.q.r = 1.0f;
    memset(oekf->state.q.i, 0, 7*sizeof(_float_t));
   // Initialize the velocity and position to 0
    memset(oekf->state.v, 0, 3*sizeof(_float_t));
    memset(oekf->state.p, 0, 3*sizeof(_float_t));
 // Initialize the covariance matrix (diagonal matrix)
 for (int i=0; i<OEKF_N; ++i) {

        for (int j=0; j<OEKF_N; ++j) {

            oekf->P[i*OEKF_N+j] = i==j ? pdiag[i] : 0;
        }
    }
}

/**
  * Runs the EKF prediction step
  * @param ekf pointer to an ekf_t structure
  * @param fx predicted values
  * @param F Jacobian of state-transition function
  * @param Q process noise matrix
  * 
  */static void ekf_predict(
        ekf_t * ekf, 
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
static void ekf_update_step3(ekf_t * ekf, _float_t GH[OEKF_N*OEKF_N])
{
    _negate(GH, OEKF_N, OEKF_N);
    _addeye(GH, OEKF_N);
    _float_t GHP[OEKF_N*OEKF_N];
    _mulmat(GH, ekf->P, GHP, OEKF_N, OEKF_N, OEKF_N);
    memcpy(ekf->P, GHP, OEKF_N*OEKF_N*sizeof(_float_t));
}

/**
  * Runs the EKF update step
  * @param ekf pointer to an ekf_t structure
  * @param z observations
  * @param hx predicted values
  * @param H sensor-function Jacobian matrix
  * @param R measurement-noise matrix
  * 
  */
static bool ekf_update(
        ekf_t * ekf, 
        const _float_t z[EKF_M], 
        const _float_t hx[OEKF_M],
        const _float_t H[EKF_M*OEKF_N],
        const _float_t R[EKF_M*EKF_M])
{        
    // G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}
    _float_t G[OEKF_N*EKF_M];
    _float_t Ht[OEKF_N*EKF_M];
    _transpose(H, Ht, EKF_M, OEKF_N);
    _float_t PHt[OEKF_N*EKF_M];
    _mulmat(ekf->P, Ht, PHt, OEKF_N, OEKF_N, EKF_M);
    _float_t HP[EKF_M*OEKF_N];
    _mulmat(H, ekf->P, HP, EKF_M, OEKF_N, OEKF_N);
    _float_t HpHt[EKF_M*EKF_M];
    _mulmat(HP, Ht, HpHt, EKF_M, OEKF_N, EKF_M);
    _float_t HpHtR[EKF_M*EKF_M];
    _addmat(HpHt, R, HpHtR, EKF_M, EKF_M);
    _float_t HPHtRinv[EKF_M*EKF_M];
    if (!invert(HpHtR, HPHtRinv)) {
        return false;
    }
    _mulmat(PHt, HPHtRinv, G, OEKF_N, EKF_M, EKF_M);

    // \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))
    _float_t z_hx[EKF_M];
    _sub(z, hx, z_hx, EKF_M);
    _float_t Gz_hx[EKF_M*OEKF_N];
    _mulvec(G, z_hx, Gz_hx, OEKF_N, EKF_M);
    _addvec(ekf->x, Gz_hx, ekf->x, OEKF_N);

    // P_k = (I - G_k H_k) P_k
    _float_t GH[OEKF_N*OEKF_N];
    _mulmat(G, H, GH, OEKF_N, EKF_M, OEKF_N);
    ekf_update_step3(ekf, GH);

    // success
    return true;
}
