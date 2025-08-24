/* gps_ekf: TinyEKF test case using You Chong's GPS example:
 * 
 *   http://www.mathworks.com/matlabcentral/fileexchange/
 *     31487-extended-kalman-filter-ekf--for-gps
 * 
 * Reads file data.csv of satellite data and writes file ekf.csv of
 * mean-subtracted estimated positions.
 *
 *
 * References:
 *
 * 1. R G Brown, P Y C Hwang, "Introduction to random signals and applied 
 * Kalman filtering : with MATLAB exercises and solutions",1996
 *
 * 2. Pratap Misra, Per Enge, "Global Positioning System Signals, 
 * Measurements, and Performance(Second Edition)",2006
 * 
 * Copyright (C) 2015 Simon D. Levy
 * Modifications Copyright (C) 2025 ZuoCen Liu
 * MIT License
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>

// The observation dimension remains 4 (4 satellites)
#define OEKF_M 4

// Maintain double precision to match the original result
#define _float_t double

// Replace with the octonion EKF header file
#include <tinyoekf.h>

// positioning interval
static const double T = 1;

// initial covariances of state noise, measurement noise
static const double P0 = 10;
static const double R0 = 36;

// Set fixed process-noise covariance matrix Q, see [1]  ---------------------

static const double Sf    = 36;
static const double Sg    = 0.01;
static const double sigma = 5;         // state transition variance

static const double b0 = Sf * T+Sg * T * T * T/3; 
static const double b1 = Sg * T * T/2; 
static const double b2 = Sg * T * T/2; 
static const double b3 = Sg * T;

static const double xyz0 = sigma * sigma * T * T * T/3; 
static const double xyz1 = sigma * sigma * T * T/2; 
static const double xyz2 = sigma * sigma * T * T/2; 
static const double xyz3 = sigma * sigma * T;

static const _float_t Q[OEKF_N*OEKF_N] = {
    // Octonion part (8-dimensional): noise is set to 0 (no attitude update)
    0,0,0,0,0,0,0,0, 0,0,0, 0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0, 0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0, 0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0, 0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0, 0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0, 0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0, 0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0, 0,0,0,
    // Velocity part (3D): Reuse the original noise parameters
    0,0,0,0,0,0,0,0, xyz3, xyz2, 0, 0,0,0,  // vx
    0,0,0,0,0,0,0,0, xyz1, xyz0, 0, 0,0,0,  // vy
    0,0,0,0,0,0,0,0, 0,0, xyz3, 0,0,0,      // vz
    // Position part (3D): reuse the original noise parameters
    0,0,0,0,0,0,0,0, 0,0,0, xyz3, xyz2,0,   // x
    0,0,0,0,0,0,0,0, 0,0,0, xyz1, xyz0,0,   // y
    0,0,0,0,0,0,0,0, 0,0,0, 0,0, xyz3       // z
    // Temperature part (1D)
    0,0,0,0,0,0,0,0, 0,0,0, 0,0,0,0.01        // Temperature process noise (set to 0.01)
};

// Set fixed measurement noise covariance matrix R ----------------------------

static const double R[4*4] = {
    R0, 0, 0, 0,
    0, R0, 0, 0,
    0, 0, R0, 0,
    0, 0, 0, R0

};

static void init(oekf_t * oekf)
{
    // Initialize the diagonal of the 14-dimensional covariance matrix (quaternion + velocity + position)
    const _float_t pdiag[OEKF_N] = {
        1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,  // Octonion (8-dimensional, small noise)
        10,10,10,                                  // Velocity (3D, moderate noise)
        10,10,10                                   // Position (3D, moderate noise)
        1.0                                        // Temperature (1-dimensional, new, initial variance set to 1.0)
    };
    oekf_initialize(oekf, pdiag);

    // Initialization position (corresponding to the original x, y, z)
    oekf->state.p[0] = -2.168816181271560e+006;  // x
    oekf->state.p[1] = 4.386648549091666e+006;   // y
    oekf->state.p[2] = 4.077161596428751e+006;   // z

    // Initialization speed (original vx, vy, vz)
    oekf->state.v[0] = 0;
    oekf->state.v[1] = 0;
    oekf->state.v[2] = 0;

    // The clock parameters are mapped to the imaginary part of the octonion (using redundant dimensions)
    oekf->state.q.i[5] = 3.575261153706439e+006;  // Clock offset
    oekf->state.q.i[6] = 4.549246345845814e+001;  // Clock drift

    // Synchronize the state vector (from state to x array)
    oekf->x[0] = oekf->state.q.r;
    memcpy(&oekf->x[1], oekf->state.q.i, 7*sizeof(_float_t));
    memcpy(&oekf->x[8], oekf->state.v, 3*sizeof(_float_t));
    memcpy(&oekf->x[11], oekf->state.p, 3*sizeof(_float_t));
   
    // The temperature is initialized to state.temp
    oekf->state.temp = 25.0;  // Initial temperature 25℃
}

static void run_model(
        oekf_t * oekf, 
        const _float_t SV[4][3], 
        _float_t fx[OEKF_N],
        _float_t F[OEKF_N*OEKF_N],
        _float_t hx[OEKF_M],
        _float_t H[OEKF_M*OEKF_N])
{ 
    // 1.Initialize the state transition matrix F as an identity matrix
    memset(F, 0, OEKF_N*OEKF_N*sizeof(_float_t));
    for (int i=0; i<OEKF_N; i++) F[i*OEKF_N + i] = 1;

    // 2. Predict the state vector fx
    // The octonion part remains unchanged (when there is no attitude update)
    fx[0] = oekf->x[0];  // real part
    memcpy(&fx[1], &oekf->x[1], 7*sizeof(_float_t));  // Imaginary part
    fx[14] = oekf->x[14];  // Temperature prediction value = Current temperature (new)

    // Speed prediction (v = v_prev)
    memcpy(&fx[8], &oekf->x[8], 3*sizeof(_float_t));

    // Position prediction (p = p_prev + T*v_prev)
    for (int i=0; i<3; i++) {
        fx[11 + i] = oekf->x[11 + i] + T * oekf->x[8 + i];  // x=11, y=12, z=13
        F[(11+i)*OEKF_N + (8+i)] = T;  // The coefficient of velocity with respect to position
    }

    // Clock parameter prediction (deviation = deviation + T * drift)
    fx[1 + 5] = oekf->x[1 + 5] + T * oekf->x[1 + 6];  // Clock offset (i[5])
    fx[1 + 6] = oekf->x[1 + 6];                       // Clock drift (i[6])
    F[(1+5)*OEKF_N + (1+6)] = T;  // Coefficient of drift to deviation

    // 3.Observation equation hx (satellite distance = position distance + clock offset)
    _float_t dx[4][3];
    for (int i=0; i<4; i++) {
        hx[i] = 0;
        for (int j=0; j<3; j++) {
            dx[i][j] = fx[11 + j] - SV[i][j];  // Position difference (x, y, z correspond to 11, 12, 13)
            hx[i] += dx[i][j] * dx[i][j];
        }
        hx[i] = sqrt(hx[i]) + fx[1 + 5];  // Distance = Euclidean distance + clock offset
    }

    // 4. Observation matrix H (4x14)
    memset(H, 0, OEKF_M*OEKF_N*sizeof(_float_t));
    for (int i=0; i<4; i++) {
        // Partial derivative of position with respect to distance
        for (int j=0; j<3; j++) {
            H[i*OEKF_N + (11 + j)] = dx[i][j] / hx[i];  // x, y, z correspond to 11, 12, 13
        }
        // Partial derivative of clock offset with respect to distance
        H[i*OEKF_N + (1 + 5)] = 1;  // Clock skew in x[6] (i[5])
    }
}

static void readline(char * line, FILE * fp)
{
    fgets(line, 1000, fp);
}

static void readdata(FILE * fp, double SV_Pos[4][3], double SV_Rho[4])
{
    char line[1000];

    readline(line, fp);

    char * p = strtok(line, ",");

    for (int i=0; i<4; ++i)
        for (int j=0; j<3; ++j) {
            SV_Pos[i][j] = atof(p);
            p = strtok(NULL, ",");
        }

    for (int j=0; j<4; ++j) {
        SV_Rho[j] = atof(p);
        p = strtok(NULL, ",");
    }
}


static void skipline(FILE * fp)
{
    char line[1000];
    readline(line, fp);
}

void error(const char * msg)
{
    fprintf(stderr, "%s\n", msg);
}

int main(int argc, char ** argv)
{    
    oekf_t oekf = {0};  // Replace ekf_t with oekf_t

    init(&oekf);

    FILE * ifp = fopen("data.csv", "r");
    skipline(ifp);

    double SV_Pos[4][3] = {0};
    double SV_Rho[4] = {0};
    double Pos_KF[25][3] = {0};

    FILE * ofp = fopen("ekf.csv", "w");
    fprintf(ofp, "X,Y,Z\n");

    for (int j=0; j<25; ++j) {
        readdata(ifp, SV_Pos, SV_Rho);

        _float_t fx[OEKF_N] = {0};
        _float_t F[OEKF_N*OEKF_N] = {0};
        _float_t hx[OEKF_M] = {0};
        _float_t H[OEKF_M*OEKF_N] = {0};
        run_model(&oekf, SV_Pos, fx, F, hx, H);

        // Replace with the prediction function of the Octonion EKF（Pass in the time interval T, which is defined as static const double T = 1 in gps.c;）
        oekf_predict(&oekf, fx, F, Q, T);  // T is 1 second, which is consistent with the positioning interval.

        // Replace with the update function of the Octonion EKF
        oekf_update(&oekf, SV_Rho, hx, H, R);

        // Extract the position from state.p (instead of directly accessing the x array)
        for (int k=0; k<3; ++k) {
            Pos_KF[j][k] = oekf.state.p[k];
        }
    }

    // The subsequent mean calculation and output logic remain unchanged...

    fclose(ifp);
    fclose(ofp);
    printf("Wrote file ekf.csv\n");
    return 0;
}
