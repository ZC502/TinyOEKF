/**
 * Custom EKF methods
 *
 * Copyright (C) 2024 Simon D. Levy
 * Modifications Copyright (C) 2025 ZuoCen Liu
 * MIT License
 */

/// @private
static void outer(
        const _float_t x[OEKF_N],
        const _float_t y[OEKF_N],
        _float_t a[OEKF_N*OEKF_N]) 
{
    for (int i=0; i<OEKF_N; i++) {
        for (int j=0; j<OEKF_N; j++) {
            a[i*OEKF_N+j] = x[i] * y[j];
        }
    }
}

/// @private
// Optimize the dot function to reduce variable accesses
static _float_t dot(const _float_t x[OEKF_N], const _float_t y[OEKF_N]) {
    _float_t d = 0;
    for (int k=0; k<OEKF_N; k++) d += x[k] * y[k];
    return d;
}

// New addition: Definition of octonion structure
typedef struct {
    _float_t r;          // real part
    _float_t i[7];       // Seven imaginary parts (i[0] to i[6])
} Octonion;

// Octonion multiplication function
static inline void octonion_mult(const Octonion *a, const Octonion *b, Octonion *c) {
    // Real part calculation
    c->r = a->r * b->r 
         - a->i[0]*b->i[0] - a->i[1]*b->i[1] - a->i[2]*b->i[2] 
         - a->i[3]*b->i[3] - a->i[4]*b->i[4] - a->i[5]*b->i[5] 
         - a->i[6]*b->i[6];

   // Imaginary part e0 (i[0]): corresponds to the standard e1, derived from the sign inversion of e2e3, e6e7, and e4e5
    c->i[0] = a->r*b->i[0] + a->i[0]*b->r 
            + (a->i[1]*b->i[2] - a->i[2]*b->i[1])   // e2e3 - e3e2 = 2e0
            + (a->i[5]*b->i[6] - a->i[6]*b->i[5])   // e6e7 - e7e6 = 2e0
            - (a->i[3]*b->i[4] - a->i[4]*b->i[3]);  // e4e5 - e5e4 = -2e0 → Subtraction is equivalent to addition

    // Imaginary part e1 (i[1]): corresponds to the standard e2, derived from the sign inversion of e3e1, e4e6, and e5e7
    c->i[1] = a->r*b->i[1] + a->i[1]*b->r 
            + (a->i[2]*b->i[0] - a->i[0]*b->i[2])   // e3e1 - e1e3 = 2e1
            + (a->i[3]*b->i[5] - a->i[5]*b->i[3])   // e4e6 - e6e4 = 2e1
            - (a->i[4]*b->i[6] - a->i[6]*b->i[4]);  // e5e7 - e7e5 = -2e1 → Subtraction is equivalent to addition

    // Imaginary part e2 (i[2]): corresponds to the standard e3, derived from the sign inversions of e1e2, e5e6, and e4e7
    c->i[2] = a->r*b->i[2] + a->i[2]*b->r 
            + (a->i[0]*b->i[1] - a->i[1]*b->i[0])   // e1e2 - e2e1 = 2e2
            + (a->i[4]*b->i[5] - a->i[5]*b->i[4])   // e5e6 - e6e5 = 2e2
            - (a->i[3]*b->i[6] - a->i[6]*b->i[3]);  // e4e7 - e7e4 = -2e2 → Subtraction is equivalent to addition

    // Imaginary part e3 (i[3]): corresponds to standard e4, derived from e1e5, e2e6, e3e7
    c->i[3] = a->r*b->i[3] + a->i[3]*b->r 
            + (a->i[0]*b->i[4] - a->i[4]*b->i[0])   // e1e5 - e5e1 = 2e3
            + (a->i[1]*b->i[5] - a->i[5]*b->i[1])   // e2e6 - e6e2 = 2e3
            + (a->i[2]*b->i[6] - a->i[6]*b->i[2]);  // e3e7 - e7e3 = 2e3

    // Imaginary part e4 (i[4]): corresponds to standard e5, derived from the sign inversion of e4e1, e2e7, e3e6
    c->i[4] = a->r*b->i[4] + a->i[4]*b->r 
            + (a->i[3]*b->i[0] - a->i[0]*b->i[3])   // e4e1 - e1e4 = 2e4
            + (a->i[1]*b->i[6] - a->i[6]*b->i[1])   // e2e7 - e7e2 = 2e4
            - (a->i[2]*b->i[5] - a->i[5]*b->i[2]);  // e3e6 - e6e3 = -2e4 → Subtraction is equivalent to addition

    // Imaginary part e5 (i[5]): corresponding to standard e6, sources e7e1, e4e2, e3e5
    c->i[5] = a->r*b->i[5] + a->i[5]*b->r 
            + (a->i[6]*b->i[0] - a->i[0]*b->i[6])   // e7e1 - e1e7 = 2e5
            + (a->i[3]*b->i[1] - a->i[1]*b->i[3])   // e4e2 - e2e4 = 2e5
            + (a->i[2]*b->i[4] - a->i[4]*b->i[2]);  // e3e5 - e5e3 = 2e5

    // Imaginary part e6 (i[6]): corresponds to standard e7, sources e1e6, e5e2, e4e3
    c->i[6] = a->r*b->i[6] + a->i[6]*b->r 
            + (a->i[0]*b->i[5] - a->i[5]*b->i[0])   // e1e6 - e6e1 = 2e6
            + (a->i[4]*b->i[1] - a->i[1]*b->i[4])   // e5e2 - e2e5 = 2e6
            + (a->i[3]*b->i[2] - a->i[2]*b->i[3]);  // e4e3 - e3e4 = 2e6
}

// Octonion normalization function
static void octonion_normalize(Octonion *q) {
    _float_t norm = sqrt(q->r*q->r + dot(q->i, q->i));  // Reuse the existing dot function
    if (norm > 0) {
        q->r /= norm;
        for (int i=0; i<7; i++) {
            q->i[i] /= norm;
        }
    }
}

// Octonion rotation vector (body frame -> navigation frame)
static void octonion_rotate(const Octonion *q, const _float_t vec_body[3], _float_t vec_nav[3]) {
    // Implementing rotation using octonion multiplication: v_nav = q * v_body * q^{-1} (simplified calculation)
    Octonion v_body_q = {0};
    memcpy(v_body_q.i, vec_body, 3*sizeof(_float_t));  // Convert a vector to a pure octonion (with a real part of 0)
    Octonion q_inv = *q;  // Octonion inverse (the inverse of a unit octonion is equal to its conjugate)
    for (int i=0; i<7; i++) q_inv.i[i] = -q_inv.i[i];
    Octonion temp, res;
    octonion_mult(q, &v_body_q, &temp);
    octonion_mult(&temp, &q_inv, &res);
    memcpy(vec_nav, res.i, 3*sizeof(_float_t));  // Take the first 3 dimensions of the imaginary part of the result
}

/**
  * Runs a custom update on the covariance matrix
  * @param ekf pointer to an oekf_t structure
  * @param A from the update P <- A P A^T
  */
static void ekf_custom_multiply_covariance(oekf_t * ekf, const _float_t A[OEKF_N*OEKF_N]) 
{
    _float_t AP[OEKF_N*OEKF_N] = {};
    _mulmat(A, ekf->P,  AP, OEKF_N, OEKF_N, OEKF_N);

    _float_t At[OEKF_N*OEKF_N] = {};
    _transpose(A, At, OEKF_N, OEKF_N);

    _mulmat(AP, At, ekf->P, OEKF_N, OEKF_N, OEKF_N);
}

/**
  * Enforces symmetry of the covariance matrix and ensures that the its values stay bounded
  * @param ekf pointer to an oekf_t structure
  * @param minval minimum covariance bound
  * @param maxval maximum covariance bound
  * 
  */
static void ekf_custom_cleanup_covariance(
        oekf_t * ekf, const float minval, const float maxval)
{

    for (int i=0; i<OEKF_N; i++) {

        for (int j=i; j<OEKF_N; j++) {

            const _float_t pval = (ekf->P[i*OEKF_N+j] + ekf->P[OEKF_N*j+i]) / 2;

            ekf->P[i*OEKF_N+j] = ekf->P[j*OEKF_N+i] =
                pval > maxval ?  maxval :
                (i==j && pval < minval) ?  minval :
                pval;
        }
    }
}

/**
  * Updates the EKF with a single scalar observation
  * @param ekf pointer to an oekf_t structure
  * @param z the observation
  * @param hx the predicted value
  * @param h one column of the sensor-function Jacobian matrix H
  * @param r one entry in the measurement-noise matrix R
  * 
  */
static void ekf_custom_scalar_update(
        oekf_t * ekf,  
        const _float_t z,
        const _float_t hx,
        const _float_t h[OEKF_N], 
        const _float_t r)

    // G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}
    _float_t ph[OEKF_N] = {};
    _mulvec(ekf->P, h, ph, OEKF_N, OEKF_N);
    const _float_t hphtr_inv = 1 / (r + dot(h, ph)); 
    _float_t g[OEKF_N] = {};
    for (int i=0; i<OEKF_N; ++i) {
        g[i] = ph[i] * hphtr_inv;
    }

    // \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))
    for (int i=0; i<OEKF_N; ++i) {
        ekf->x[i] += g[i] * (z - hx);
    }

    // P_k = (I - G_k H_k) P_k$
    _float_t GH[OEKF_N*OEKF_N];
    outer(g, h, GH); 
    ekf_update_step3(ekf, GH);

    // Does this belong here, or in caller?
    for (int i=0; i<OEKF_N; i++) {
        for (int j=i; j<OEKF_N; j++) {
            ekf->P[i*OEKF_N+j] += r * g[i] * g[j];
        }
    }
}
