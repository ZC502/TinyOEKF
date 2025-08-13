/**
 * Custom EKF methods
 *
 * Copyright (C) 2024 Simon D. Levy
 * Modifications Copyright (C) 2025 ZuoCen Liu
 * MIT License
 */

/// @private
static void outer(
        const _float_t x[EKF_N],
        const _float_t y[EKF_N],
        _float_t a[EKF_N*EKF_N]) 
{
    for (int i=0; i<EKF_N; i++) {
        for (int j=0; j<EKF_N; j++) {
            a[i*EKF_N+j] = x[i] * y[j];
        }
    }
}

/// @private
static _float_t dot(const _float_t x[EKF_N], const _float_t y[EKF_N]) 
{
    _float_t d = 0;

    for (int k=0; k<EKF_N; k++) {
        d += x[k] * y[k];
    }

    return d;
}

// New addition: Definition of octonion structure
typedef struct {
    _float_t r;          // real part
    _float_t i[7];       // Seven imaginary parts (i[0] to i[6])
} Octonion;

// New: Octonion multiplication function
static inline void octonion_mult(const Octonion *a, const Octonion *b, Octonion *c) {
    // Real part calculation
    c->r = a->r * b->r 
         - a->i[0]*b->i[0] - a->i[1]*b->i[1] - a->i[2]*b->i[2] 
         - a->i[3]*b->i[3] - a->i[4]*b->i[4] - a->i[5]*b->i[5] 
         - a->i[6]*b->i[6];

    // Imaginary part e0 (i[0])
    c->i[0] = a->r*b->i[0] + a->i[0]*b->r 
            + (a->i[1]*b->i[3] - a->i[3]*b->i[1]) 
            + (a->i[2]*b->i[6] - a->i[6]*b->i[2]) 
            + (a->i[4]*b->i[2] - a->i[2]*b->i[4]) 
            + (a->i[5]*b->i[1] - a->i[1]*b->i[5]) 
            + (a->i[6]*b->i[5] - a->i[5]*b->i[6]);

     // Imaginary part e1 (i[1])
    c->i[1] = a->r*b->i[1] + a->i[1]*b->r 
            + (a->i[0]*b->i[3] - a->i[3]*b->i[0]) 
            + (a->i[2]*b->i[4] - a->i[4]*b->i[2]) 
            + (a->i[5]*b->i[2] - a->i[2]*b->i[5]) 
            + (a->i[6]*b->i[5] - a->i[5]*b->i[6]) 
            + (a->i[0]*b->i[5] - a->i[5]*b->i[0]);

    // Imaginary part e2 (i[2])
    c->i[2] = a->r*b->i[2] + a->i[2]*b->r 
            + (a->i[0]*b->i[6] - a->i[6]*b->i[0]) 
            + (a->i[1]*b->i[4] - a->i[4]*b->i[1]) 
            + (a->i[3]*b->i[5] - a->i[5]*b->i[3]) 
            + (a->i[4]*b->i[0] - a->i[0]*b->i[4]) 
            + (a->i[5]*b->i[1] - a->i[1]*b->i[5]);

    // Imaginary part e3 (i[3])
    c->i[3] = a->r*b->i[3] + a->i[3]*b->r 
            + (a->i[1]*b->i[0] - a->i[0]*b->i[1]) 
            + (a->i[2]*b->i[5] - a->i[5]*b->i[2]) 
            + (a->i[4]*b->i[1] - a->i[1]*b->i[4]) 
            + (a->i[6]*b->i[4] - a->i[4]*b->i[6]) 
            + (a->i[0]*b->i[6] - a->i[6]*b->i[0]);

    // Imaginary part e4 (i[4])
    c->i[4] = a->r*b->i[4] + a->i[4]*b->r 
            + (a->i[1]*b->i[2] - a->i[2]*b->i[1]) 
            + (a->i[3]*b->i[1] - a->i[1]*b->i[3]) 
            + (a->i[4]*b->i[5] - a->i[5]*b->i[4]) 
            + (a->i[6]*b->i[3] - a->i[3]*b->i[6]) 
            + (a->i[0]*b->i[2] - a->i[2]*b->i[0]);

    // Imaginary part e5 (i[5])
    c->i[5] = a->r*b->i[5] + a->i[5]*b->r 
            + (a->i[1]*b->i[2] - a->i[2]*b->i[1]) 
            + (a->i[2]*b->i[3] - a->i[3]*b->i[2]) 
            + (a->i[4]*b->i[6] - a->i[6]*b->i[4]) 
            + (a->i[6]*b->i[0] - a->i[0]*b->i[6]) 
            + (a->i[0]*b->i[1] - a->i[1]*b->i[0]);

    // Imaginary part e6 (i[6])
    c->i[6] = a->r*b->i[6] + a->i[6]*b->r 
            + (a->i[0]*b->i[2] - a->i[2]*b->i[0]) 
            + (a->i[1]*b->i[5] - a->i[5]*b->i[1]) 
            + (a->i[4]*b->i[5] - a->i[5]*b->i[4]) 
            + (a->i[3]*b->i[6] - a->i[6]*b->i[3]) 
            + (a->i[0]*b->i[3] - a->i[3]*b->i[0]);
}

/**
  * Runs a custom update on the covariance matrix
  * @param ekf pointer to an ekf_t structure
  * @param A from the update P <- A P A^T
  */
static void ekf_custom_multiply_covariance(
        ekf_t * ekf, const _float_t A[EKF_N*EKF_N]) 
{
    _float_t AP[EKF_N*EKF_N] = {};
    _mulmat(A, ekf->P,  AP, EKF_N, EKF_N, EKF_N);

    _float_t At[EKF_N*EKF_N] = {};
    _transpose(A, At, EKF_N, EKF_N);

    _mulmat(AP, At, ekf->P, EKF_N, EKF_N, EKF_N);
}

/**
  * Enforces symmetry of the covariance matrix and ensures that the its values stay bounded
  * @param ekf pointer to an ekf_t structure
  * @param minval minimum covariance bound
  * @param maxval maximum covariance bound
  * 
  */
static void ekf_custom_cleanup_covariance(
        ekf_t * ekf, const float minval, const float maxval)
{

    for (int i=0; i<EKF_N; i++) {

        for (int j=i; j<EKF_N; j++) {

            const _float_t pval = (ekf->P[i*EKF_N+j] + ekf->P[EKF_N*j+i]) / 2;

            ekf->P[i*EKF_N+j] = ekf->P[j*EKF_N+i] =
                pval > maxval ?  maxval :
                (i==j && pval < minval) ?  minval :
                pval;
        }
    }
}

/**
  * Updates the EKF with a single scalar observation
  * @param ekf pointer to an ekf_t structure
  * @param z the observation
  * @param hx the predicted value
  * @param h one column of the sensor-function Jacobian matrix H
  * @param r one entry in the measurement-noise matrix R
  * 
  */
static void ekf_custom_scalar_update(
        ekf_t * ekf,
        const _float_t z,
        const _float_t hx,
        const _float_t h[EKF_N], 
        const _float_t r)
{
    (void)ekf_update;

    // G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}
    _float_t ph[EKF_N] = {};
    _mulvec(ekf->P, h, ph, EKF_N, EKF_N);
    const _float_t hphtr_inv = 1 / (r + dot(h, ph)); 
    _float_t g[EKF_N] = {};
    for (int i=0; i<EKF_N; ++i) {
        g[i] = ph[i] * hphtr_inv;
    }

    // \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))
    for (int i=0; i<EKF_N; ++i) {
        ekf->x[i] += g[i] * (z - hx);
    }

    // P_k = (I - G_k H_k) P_k$
    _float_t GH[EKF_N*EKF_N];
    outer(g, h, GH); 
    ekf_update_step3(ekf, GH);

    // Does this belong here, or in caller?
    for (int i=0; i<EKF_N; i++) {
        for (int j=i; j<EKF_N; j++) {
            ekf->P[i*EKF_N+j] += r * g[i] * g[j];
        }
    }
}
