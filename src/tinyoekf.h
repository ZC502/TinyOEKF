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
#include "tinyoekf_custom.h"  // Include the header file of custom functions

#ifndef _float_t
#define _float_t float
#endif

#define OEKF_N 14  // 8(octonions attitude) + 3(velocity) + 3(position)
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
static void oekf_predict(oekf_t *ekf, const _float_t fx[OEKF_N], const _float_t F[OEKF_N * OEKF_N], const _float_t Q[OEKF_N * OEKF_N],
const _float_t dt, const _float_t accel_body[3])  // Add the dt parameter,Airframe acceleration parameters (sensor input)
{

// Default octonion motion model (when fx or F is NULL)
    _float_t default_fx[OEKF_N];
    _float_t default_F[OEKF_N*OEKF_N] = {0};
    if (fx == NULL || F == NULL) {

        // 1. Obtain angular velocity from the sensor (not the state vector, as the state does not contain angular velocity)
        _float_t omega[3] = {...};   // Passed in from the outside or received through parameters

        // 2. Octonion update (nonlinear attitude modeling based on angular velocity)
        Octonion dq;  // Incremental octonion
               _float_t theta = sqrt(omega[0]*omega[0] + omega[1]*omega[1] + omega[2]*omega[2]) * dt;
        dq.r = cos(theta/2);
        dq.i[0] = omega[0] * sin(theta/2) / (theta + 1e-8);  // Avoid division by zero
        dq.i[1] = omega[1] * sin(theta/2) / (theta + 1e-8);
        dq.i[2] = omega[2] * sin(theta/2) / (theta + 1e-8);
        memset(&dq.i[3], 0, 4*sizeof(_float_t));  // The imaginary part of high dimensions is not used for the time being.
        octonion_mult(&ekf->state.q, &dq, &q_new);  // Nonlinear attitude update

        // 3. Velocity update (coupled attitude: conversion of body acceleration to navigation frame)
        _float_t accel_nav[3];  // Navigation system acceleration
        octonion_rotate(&ekf->state.q, accel_body, accel_nav);  // Add a new attitude rotation function
        for (int i=0; i<3; i++) {
            default_fx[8+i] = ekf->x[8+i] + (accel_nav[i] - 9.81) * dt;  // Subtract gravity
        }
        
        // 4. Position Update (Coupled Velocity)
        for (int i=0; i<3; i++) {
            default_fx[11+i] = ekf->x[11+i] + default_fx[8+i] * dt;  // Position = Position + Velocity * dt
        }
        // Update octonion: q_new = q_old * dq (non-commutative multiplication)
        Octonion q_new;
        octonion_mult(&ekf->state.q, &dq, &q_new);
        octonion_normalize(&q_new);
        // Mapped to default_fx
        default_fx[0] = q_new.r;
        memcpy(&default_fx[1], q_new.i, 7*sizeof(_float_t));
        // Velocity/position update (v = v + adt, p = p + vdt)
        for (int i=0; i<3; i++) {
            default_fx[8+i] = ekf->x[8+i] + ekf->x[17+i] * dt;  // Assume that x[17-19] is acceleration
            default_fx[11+i] = ekf->x[11+i] + default_fx[8+i] * dt;
        }
        // Generate the default F matrix (simplified version, with only the diagonal elements being 1)
        for (int i=0; i<OEKF_N; i++) default_F[i*OEKF_N +i] = 1;
        fx = default_fx;
        F = default_F;
    }

// Update the state vector
    memcpy(ekf->x, fx, OEKF_N * sizeof(_float_t));
    
    // Update the state structure (synchronized with vector x)
    ekf->state.q.r = ekf->x[0];
    memcpy(ekf->state.q.i, &ekf->x[1], 7 * sizeof(_float_t));
    memcpy(ekf->state.v, &ekf->x[8], 3 * sizeof(_float_t));
    memcpy(ekf->state.p, &ekf->x[11], 3 * sizeof(_float_t));
    
    // Normalize the octonion to ensure the validity of the state
    octonion_normalize(&ekf->state.q);
    
    // Calculate the covariance matrix P = F * P * F^T + Q
    _float_t FP[OEKF_N * OEKF_N] = {};
    _mulmat(F, ekf->P, FP, OEKF_N, OEKF_N, OEKF_N);
    
    _float_t Ft[OEKF_N * OEKF_N] = {};
    _transpose(F, Ft, OEKF_N, OEKF_N);
    
    _float_t FPFt[OEKF_N * OEKF_N] = {};
    _mulmat(FP, Ft, FPFt, OEKF_N, OEKF_N, OEKF_N);
    
    _addmat(FPFt, Q, ekf->P, OEKF_N, OEKF_N);
}

// Residual Detection and Noise Adjustment Function
static bool oekf_detect_perturbation(const _float_t z[OEKF_M], const _float_t hx[OEKF_M], const _float_t threshold) {
    _float_t res[OEKF_M];
    _sub(z, hx, res, OEKF_M); // Reuse the existing _sub function to calculate residuals
    _float_t res_norm = 0;
    for (int i=0; i<OEKF_M; i++) res_norm += res[i]*res[i];
    return res_norm > threshold*threshold;
}

static void oekf_adapt_Q(_float_t Q[OEKF_N*OEKF_N], const bool is_perturbed, const _float_t scale) {
    if (is_perturbed) {
        for (int i=0; i<OEKF_N*OEKF_N; i++) Q[i] *= scale;
    }
}

// Asynchronous prediction function with timestamp alignment
static void oekf_predict_async(oekf_t *ekf, const _float_t dt, const _float_t Q[OEKF_N*OEKF_N]) {
    // Call the predict function with dt, using the default motion model
    oekf_predict(ekf, NULL, NULL, Q, dt);
}

// Update step helper function
static void oekf_update_step3(oekf_t *ekf, _float_t GH[OEKF_N * OEKF_N]) {
    _negate(GH, OEKF_N, OEKF_N);
    _addeye(GH, OEKF_N);
    _float_t GHP[OEKF_N * OEKF_N];
    _mulmat(GH, ekf->P, GHP, OEKF_N, OEKF_N, OEKF_N);
    memcpy(ekf->P, GHP, OEKF_N * OEKF_N * sizeof(_float_t));
}

/**
 * Lightweight update for small-dimensional observations (m ≤ 3) to avoid complete matrix inversion
 * @param Pointer to the ekf oekf_t structure
 * @param z Observation vector (length m)
 * @param hx Predicted observed value (length m)
 * @param H Observation matrix (m x OEKF_N)
 * @param R Observation noise matrix (m x m, diagonal matrix)
 * @param m Observation dimension (such as 3D IMU)
 */
static bool oekf_update_small(
        oekf_t *ekf, 
        const _float_t z[], 
        const _float_t hx[], 
        const _float_t H[], 
        const _float_t R[], 
        const int m) 
{
    // 1. Calculate the residual y = z - hx
    _float_t y[m];
    _sub(z, hx, y, m);

    // 2. Calculate P*H^T (OEKF_N × m)
    _float_t PHt[OEKF_N * m];
    _mulmat(ekf->P, H, PHt, OEKF_N, OEKF_N, m);  // H is stored column-wise at this point, which is equivalent to H^T

    // 3. Calculate HPH^T + R (m x m, simplified using the properties of diagonal matrices)
    _float_t S[m * m];
    memset(S, 0, m*m*sizeof(_float_t));
    for (int i = 0; i < m; i++) {
        // Diagonal elements: Dot product of row i of H and column i of PHt + R[i][i]
        S[i*m + i] = dot(&H[i*OEKF_N], &PHt[i], OEKF_N) + R[i*m + i];
        // Off-diagonal elements are 0 (assuming R is a diagonal matrix, commonly used for small-dimensional observations)
    }

    // 4. Calculate the Kalman gain G = P*H^T * S^{-1} (using the diagonal property of S, directly take the reciprocal)
    _float_t G[OEKF_N * m];
    for (int i = 0; i < OEKF_N; i++) {
        for (int j = 0; j < m; j++) {
            G[i*m + j] = PHt[i*m + j] / S[j*m + j];  // The inverse of S is the reciprocal of the diagonal elements.
        }
    }

    // 5. Update state x = x + G*y
    _float_t Gy[OEKF_N];
    _mulvec(G, y, Gy, OEKF_N, m);
    _addvec(ekf->x, Gy, ekf->x, OEKF_N);

    // 6. Synchronize the state to the structure and normalize the quaternion
    ekf->state.q.r = ekf->x[0];
    memcpy(ekf->state.q.i, &ekf->x[1], 7*sizeof(_float_t));
    memcpy(ekf->state.v, &ekf->x[8], 3*sizeof(_float_t));
    memcpy(ekf->state.p, &ekf->x[11], 3*sizeof(_float_t));
    octonion_normalize(&ekf->state.q);

    // 7. Update covariance P = (I - G*H) * P (lightweight implementation)
    _float_t GH[OEKF_N * OEKF_N];
    _mulmat(G, H, GH, OEKF_N, m, OEKF_N);
    oekf_update_step3(ekf, GH);  // Reuse the existing covariance update logic
    return true;
}

// Update steps
// Modify the function parameters, add the observation dimension m, and remove the fixed OEKF_M
static bool oekf_update(
        oekf_t *ekf, 
        const _float_t z[], 
        const _float_t hx[], 
        const _float_t H[], 
        const _float_t R[], 
        const int m)  // Observation dimension (e.g., 3 for IMU)
{
 // 对小维度观测（m <= 3）使用轻量化更新，跳过矩阵求逆
    if (m > 0 && m <= 3) {
        return oekf_update_small(ekf, z, hx, H, R, m);
    }
    // Calculate the gain G = P * H^T * (H * P * H^T + R)^-1
    _float_t G[OEKF_N * m];
    _float_t Ht[OEKF_N * m;
    _transpose(H, Ht, m, OEKF_N);
    
    _float_t PHt[OEKF_N * m];
    _mulmat(ekf->P, Ht, PHt, OEKF_N, OEKF_N, OEKF_M);
    
    _float_t HP[OEKF_M * OEKF_N];
    _mulmat(H, ekf->P, HP, m, OEKF_N, OEKF_N);
    
    _float_t HpHt[OEKF_M * OEKF_M];
    _mulmat(HP, Ht, HpHt, m, OEKF_N, m);
    
    _float_t HpHtR[OEKF_M * OEKF_M];
    _addmat(HpHt, R, HpHtR, m, m);
    
    _float_t HPHtRinv[OEKF_M * OEKF_M];
    if (!invert(HpHtR, HPHtRinv)) {
        return false;
    }
    
    _mulmat(PHt, HPHtRinv, G, OEKF_N, m, m);
    
    // Update the state vector x = x + G*(z - hx)
    _float_t z_hx[m];
    _sub(z, hx, z_hx, m);

    // Added: Detect disturbances and adjust Q (if adjustment is needed before updating)
    _float_t adapted_Q[OEKF_N*OEKF_N];
    
    _float_t Gz_hx[OEKF_N];
    _mulvec(G, z_hx, Gz_hx, OEKF_N, m);
    
    _addvec(ekf->x, Gz_hx, ekf->x, OEKF_N);
    
    // Update the state structure (synchronized with vector x)
    ekf->state.q.r = ekf->x[0];
    memcpy(ekf->state.q.i, &ekf->x[1], 7 * sizeof(_float_t));
    memcpy(ekf->state.v, &ekf->x[8], 3 * sizeof(_float_t));
    memcpy(ekf->state.p, &ekf->x[11], 3 * sizeof(_float_t));
    // Normalize the octonion to ensure the validity of the state
    octonion_normalize(&ekf->state.q);

    // Geometric constraints (positional/velocity physical limitations)
    // Position constraint (assuming the z-axis is the height, with a minimum value limit of 0)
    if (ekf->state.p[2] < 0) ekf->state.p[2] = 0;
    // Speed constraint (limiting the maximum speed, such as 20m/s)
    _float_t v_max = 20;
    for (int i=0; i<3; i++) {
        if (ekf->state.v[i] > v_max) ekf->state.v[i] = v_max;
        else if (ekf->state.v[i] < -v_max) ekf->state.v[i] = -v_max;
    }
    // Synchronize the constrained state to the x vector
    memcpy(&ekf->x[11], ekf->state.p, 3*sizeof(_float_t));
    memcpy(&ekf->x[8], ekf->state.v, 3*sizeof(_float_t));
    
    // Update the covariance matrix P = (I - G*H) * P
    _float_t GH[OEKF_N * OEKF_N];
    _mulmat(G, H, GH, OEKF_N, m, OEKF_N);
    oekf_update_step3(ekf, GH);
    
    return true;
}

#endif // TINYOEKF_H

/// @private
static void _mulmat(
        const _float_t * a, 
        const _float_t * b, 
        _float_t * c, 
        const int arows, 
        const int acols, 
        const int bcols)
{
    memset(c, 0, arows*bcols*sizeof(_float_t));
    
    // Optimization is only performed for the expansion of 14x14 matrices (to adapt to the scenario where OEKF_N=14)
    if (arows == 14 && acols == 14 && bcols == 14) {
        // Line 0 expansion：c[0][j] = sum(k=0~13) a[0][k] * b[k][j]
        for (int j = 0; j < 14; j++) {
            c[0*14 + j] = 
                a[0*14 + 0] * b[0*14 + j] +
                a[0*14 + 1] * b[1*14 + j] +
                a[0*14 + 2] * b[2*14 + j] +
                a[0*14 + 3] * b[3*14 + j] +
                a[0*14 + 4] * b[4*14 + j] +
                a[0*14 + 5] * b[5*14 + j] +
                a[0*14 + 6] * b[6*14 + j] +
                a[0*14 + 7] * b[7*14 + j] +
                a[0*14 + 8] * b[8*14 + j] +
                a[0*14 + 9] * b[9*14 + j] +
                a[0*14 + 10] * b[10*14 + j] +
                a[0*14 + 11] * b[11*14 + j] +
                a[0*14 + 12] * b[12*14 + j] +
                a[0*14 + 13] * b[13*14 + j];
        }

        // Line 1 expansion:c[1][j] = sum(k=0~13) a[1][k] * b[k][j]
        for (int j = 0; j < 14; j++) {
            c[1*14 + j] = 
                a[1*14 + 0] * b[0*14 + j] +
                a[1*14 + 1] * b[1*14 + j] +
                a[1*14 + 2] * b[2*14 + j] +
                a[1*14 + 3] * b[3*14 + j] +
                a[1*14 + 4] * b[4*14 + j] +
                a[1*14 + 5] * b[5*14 + j] +
                a[1*14 + 6] * b[6*14 + j] +
                a[1*14 + 7] * b[7*14 + j] +
                a[1*14 + 8] * b[8*14 + j] +
                a[1*14 + 9] * b[9*14 + j] +
                a[1*14 + 10] * b[10*14 + j] +
                a[1*14 + 11] * b[11*14 + j] +
                a[1*14 + 12] * b[12*14 + j] +
                a[1*14 + 13] * b[13*14 + j];
        }

        // Keep the loop in lines 2 to 13 (to balance code amount and performance)
        for (int i = 2; i < 14; i++) {
            for (int j = 0; j < 14; j++) {
                for (int k = 0; k < 14; k++) {
                    c[i*14 + j] += a[i*14 + k] * b[k*14 + j];
                }
            }
        }
    }
    // Non-14x14 matrices use regular loops (compatible with other dimensions)
    else {
        for (int i = 0; i < arows; i++) {
            for (int j = 0; j < bcols; j++) {
                for (int k = 0; k < acols; k++) {
                    c[i*bcols + j] += a[i*acols + k] * b[k*bcols + j];
                }
            }
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
 * @param ekf pointer to an oekf_t structure
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
  * @param ekf pointer to an oekf_t structure
  * @param fx predicted values
  * @param F Jacobian of state-transition function
  * @param Q process noise matrix
  * 
  */static void ekf_predict(
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
static void ekf_update_step3(oekf_t * ekf, _float_t GH[OEKF_N*OEKF_N])
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
static bool ekf_update(
        oekf_t * ekf, 
        const _float_t z[OEKF_M], 
        const _float_t hx[OEKF_M],
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
    _addmat(HpHt, R, HpHtR, OEKF_M,OEKF_M);
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
    ekf_update_step3(ekf, GH);

    // success
    return true;
}
