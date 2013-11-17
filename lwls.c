#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A11 d0
#define A12 d1
#define A21 d1
#define A13 d2
#define A22 d2
#define A31 d2
#define A14 d3
#define A23 d3
#define A32 d3
#define A41 d3
#define A24 d4
#define A33 d4
#define A42 d4
#define A34 d5
#define A43 d5
#define A44 d6

// typedef enum {
//   epan,
//   rect,
//   gauss,
//   gausvar,
//   quar
// } kernel_type;

// int bin_subjects (double *x, double *y, int *num_bins,
//   int *n, int *ni, double *ans) {
//   ns = *n;
//   for ()
// }

int lwls(
         // input stuff
         double *bandwidthP, int *kernelP, double *x_in, double *y_in, double *count_in,
         // data describing stuff
         int *n_inP, int *n_outP, int *weight_out_need,
         // output stuff
         double *x_out, double *output, double *weight_out,
         // derivative and cross validation stuff
         int *powerP, int *cv_modeP, double *cv_value) {
    
    // power: 1 mean
    //        2 1st derivative
    //        3 2nd derivative
    
    double bandwidth = *bandwidthP;
    int kernel = *kernelP, n_in = *n_inP, n_out = *n_outP, power = *powerP, cv_mode = *cv_modeP;
    
    // cv_mode: 0 non-cv
    //          1 ocv
    //          2 gcv
    //          3 geometric mean
    // When cv mode is on, n_out should equal to n_in
    
    // declaration for cv
    double *h = NULL, *h_init = NULL; // h_{ii} or trH / n
    double *mu = NULL, *mu_init = NULL;  // for mean function used in cv, maybe a better name?
    double sumh = 0;
    if (cv_mode != 0) {
        // when do cv, we need save mu somewhere
        // when power is 1, mu is same to output
        // but in other case, we need some place to save mu separately
        // power = 1; //set output to mu
        h = (double *) calloc(n_out, sizeof(double));
        h_init = h;
        mu = (double *) calloc(n_out, sizeof(double));
        mu_init = mu;
    }
    // end of declaration for cv
    
    // declaration for 2d weight
    double *temp_weight = NULL, *temp_weight_init = NULL; // save weight for 1 point temporarily
    if (*weight_out_need) {
        temp_weight = (double *) calloc(n_in, sizeof(double));
        temp_weight_init = temp_weight;
    }
    double *weight_outP = NULL;
    // end of declaration for 2d weight
    
    double d0 = 0, d1 = 0, d2 = 0, d3 = 0, d4 = 0, d5 = 0, d6 = 0,
    b1 = 0, b2 = 0, b3 = 0, b4 = 0, cur_x_out = 0, d = 0, d_bw = 0, w = 0, min_bw = x_in[n_in - 1];
    int in_index = 0, out_index = 0, i = 0;
    
    if (bandwidth < 0) {
        fprintf(stderr, "Bandwidth choice for mu(t) and/or its derivative must be positive!\n");
        return 1;
    }
    // loop over out point, set output and weight for next lwls seq
    for (out_index = 0; out_index < n_out; out_index++) {
        d0 = d1 = d2 = d3 = d4 = d5 = d6 = b1 = b2 = b3 = b4 = 0;
        cur_x_out = x_out[out_index];
        for (in_index = 0; in_index < n_in; in_index++) {
            // May have better check choice, like isnan(y_in[in_index])
            if (count_in[in_index] == 0) {
                continue;
            }
            if (kernel == 0 || kernel == 1 || kernel == 4) {
                // when using epan, rect, or quar
                // using points within 1 bandwidth
                if (x_in[in_index] < cur_x_out - bandwidth) {
                    continue;
                }
                if (x_in[in_index] > cur_x_out + bandwidth) {
                    // reweight
                    if (*weight_out_need) {
                        weight_outP = &(weight_out[out_index]);
                        for (i = 0; i != n_in; i++) {
                            (*weight_outP) += temp_weight[i] * temp_weight[i] / count_in[i];
                            temp_weight[i] = 0;
                        }
                        (*weight_outP) = 1 / (*weight_outP);
                    }
                    break;
                }
            }
            else {
                // when using gauss, gausvar,
                // using points within 5 bandwidth
                // accurate enough
                if (x_in[in_index] < cur_x_out - 5 * bandwidth) {
                    continue;
                }
                if (x_in[in_index] > cur_x_out + 5 * bandwidth) {
                    // reweight
                    if (*weight_out_need) {
                        weight_outP = &(weight_out[out_index]);
                        for (i = 0; i != n_in; i++) {
                            (*weight_outP) += temp_weight[i] * temp_weight[i] / count_in[i];
                            temp_weight[i] = 0;
                        }
                        (*weight_outP) = 1 / (*weight_outP);
                    }
                    break;
                }
            }
            d = cur_x_out - x_in[in_index];
            d_bw = d / bandwidth;
            switch (kernel) {
                    // here w is the real weight
                case 0:
                    w = count_in[in_index] * (1 - d_bw * d_bw);
                    break;
                case 1:
                    w = count_in[in_index];
                    break;
                case 2:
                    w = count_in[in_index] / exp(d_bw * d_bw / 2);
                    break;
                case 3:
                    d_bw = d_bw * d_bw;
                    w = count_in[in_index] / exp(d_bw / 2) * (1.25 - 0.25 * d_bw);
                    break;
                case 4:
                    w = count_in[in_index] * (1 - 2 * d_bw * d_bw + d_bw * d_bw * d_bw * d_bw);
                    break;
            }
            // prepare weight for 2d smoothing
            if (*weight_out_need) {
                temp_weight[in_index] = w;
            }
            d0 += w;
            d1 += w * d;
            d2 += w * d * d;
            b1 += w * y_in[in_index];
            b2 += w * d * y_in[in_index];
            if (power > 1) {
                d3 += w * d * d * d;
                d4 += w * d * d * d * d;
                b3 += w * d * d * y_in[in_index];
            }
            if (power > 2) {
                d5 += w * d * d * d * d * d;
                d6 += w * d * d * d * d * d * d;
                b4 += w * d * d * d * y_in[in_index];
            }
            // reweight
            if (*weight_out_need && in_index == n_in - 1) {
                weight_outP = &(weight_out[out_index]);
                for (i = 0; i != n_in; i++) {
                    (*weight_outP) += temp_weight[i] * temp_weight[i] / count_in[i];
                    temp_weight[i] = 0;
                }
                (*weight_outP) = 1 / (*weight_outP);
            }
        }
        // calculating cv value
        if (power == 1) {
            output[out_index] = (b1*d2 - b2*d1)/(- d1*d1 + d0*d2);
            if (cv_mode) {
                //set h_{ii} and mu
                h[out_index] = count_in[out_index] * A22 / (A11 * A22 - A12 * A12);
                mu[out_index] = output[out_index];
            }
        }
        if (power == 2) {
            output[out_index] = -(b2*d2*d2 + b1*d1*d4 - b1*d2*d3 - b2*d0*d4 + b3*d0*d3 - b3*d1*d2)/(d4*d1*d1 - 2*d1*d2*d3 + d2*d2*d2 - d0*d4*d2 + d0*d3*d3);
            if (cv_mode) {
                //set h_{ii} and mu
                h[out_index] = count_in[out_index] * A22 / (A11 * A22 - A12 * A12);
                //**************************************************//
                mu[out_index] = 0;
            }
        }
        if (power == 3) {
            output[out_index] = 2*(b3*d0*d4*d4 - b1*d2*d4*d4 - b2*d3*d3*d3 + b1*d3*d3*d4 + b3*d2*d3*d3 + b4*d1*d3*d3 + b1*d2*d2*d6 - b4*d2*d2*d3 + b3*d1*d1*d6 - b4*d1*d1*d5 - b1*d1*d3*d6 + b1*d1*d4*d5 - b1*d2*d3*d5 + b2*d0*d3*d6 - b2*d0*d4*d5 - b2*d1*d2*d6 + b2*d1*d3*d5 + b2*d2*d3*d4 - b3*d0*d2*d6 - 2*b3*d1*d3*d4 + b4*d0*d2*d5 - b4*d0*d3*d4 + b4*d1*d2*d4)/(d6*d1*d1*d4 - d1*d1*d5*d5 - 2*d6*d1*d2*d3 + 2*d1*d2*d4*d5 + 2*d1*d3*d3*d5 - 2*d1*d3*d4*d4 + d6*d2*d2*d2 - 2*d2*d2*d3*d5 - d2*d2*d4*d4 + 3*d2*d3*d3*d4 - d0*d6*d2*d4 + d0*d2*d5*d5 - d3*d3*d3*d3 + d0*d6*d3*d3 - 2*d0*d3*d4*d5 + d0*d4*d4*d4);
            if (cv_mode) {
                //set h_{ii} and mu
                h[out_index] = count_in[out_index] * A22 / (A11 * A22 - A12 * A12);
                //**************************************************//
                mu[out_index] = 0;
            }
        }
    }
    //end of out_index loop
    
    // gcv or geomean
    if (cv_mode > 1) {
        for (i = 0; i < n_in; i++) {
            sumh += h[i];
        }
        sumh /= n_in;
        for (i = 0; i < n_in; i++) {
            h[i] = sumh;
        }
    }
    
    // assign cv value
    if (cv_mode) {
        *cv_value = 0;
        for (i = 0; i < n_in; i++) {
            (*cv_value) += (y_in[i] - mu[i]) * (y_in[i] - mu[i]) / (1 - h[i]) / (1 - h[i]);
        }
        (*cv_value) /= n_in;
        if (isnan(*cv_value)) {
            *cv_value = INFINITY;
        }
    }
    
    // clean memory
    if (h != NULL) {
        free(h_init);
        free(mu_init);
        h = NULL;
        mu = NULL;
        h_init = NULL;
        mu_init = NULL;
    }
    if (temp_weight != NULL) {
        free(temp_weight_init);
        temp_weight = NULL;
        temp_weight_init = NULL;
    }
    return 0;
}

//int binning2d(double *x1, double *x2, double *y, int *num_raw,
//              double *a10, double *a20, double *bin_length1, double *bin_length2,
//              int *num_bins1, int *num_bins2, double *output, int *count) {
//    int nb1 = num_bins1[0], ind = 0;
//    int row_ind = 0, col_ind = 0;
//    double a1s = a10[0], a2s = a20[0], bl1 = bin_length1[0], bl2 = bin_length2[0];
//    for (int i = 0; i != *num_raw; i++) {
//        row_ind = (int) ((x2[i] - a2s) / bl2);
//        col_ind = (int) ((x1[i] - a1s) / bl1);
//        if (row_ind > (num_bins2[0] - 1)) {
//            row_ind--;
//        }
//        if (col_ind > (num_bins1[0] - 1)) {
//            col_ind--;
//        }
//        ind = row_ind * nb1 + col_ind;
//        output[ind] += y[i];
//        (count[ind])++;
//    }
//    for (int i = 0; i != num_bins1[0] * num_bins2[0]; i++) {
//        if (count[i] != 0) {
//            output[i] /= count[i];
//        }
//    }
//    return 0;
//}


int lwls_seq (
              // input stuff
              // bandwidth and kernel are same for all columns
              // x_in: vector with length n
              // y_in: vector with length n * m
              // w_in: vector with length n * m
              //       in 2d smoothing, use count_in in 1st step
              //       and weight_out in 2nd step
              double *bandwidthP, int *kernelP, double *x_in, double *y_in, double *w_in, 
              // data describing stuff
              // n_inP:  number of y
              // n_outP: number of output
              // n_colP: number of columns
              int *n_inP, int *n_outP, int *n_colP, int *weight_out_need, 
              // output stuff
              // x_out:  vector with length n'
              // output: vector with length n' * m
              double *x_out, double *output, double *weight_out, 
              // derivative and cross validation stuff
              int *powerP, int *cv_modeP, double *cv_value) {
    int n_row_in = *n_inP / *n_colP, n_row_out = *n_outP / *n_colP;
    if (n_row_in * (*n_colP) != *n_inP) {
        fprintf(stderr, "Number of input is not a multiple of number of column!\n");
        return 1;
    }
    if (n_row_out * (*n_colP) != *n_inP) {
        fprintf(stderr, "Number of output is not a multiple of number of column!\n");
        return 1;
    }
    int col_index = 0;
    double *y_col_in = y_in, *w_col_in = w_in, *output_col = output, *w_col_out = weight_out;
    double sum_cv_value = 0;
    for (col_index = 0; col_index != *n_colP; col_index++) {
        lwls(bandwidthP, kernelP, x_in, y_col_in, w_col_in,
             &n_row_in, &n_row_out, weight_out_need, 
             x_out, output_col, w_col_out,
             powerP, cv_modeP, cv_value);
        printf("%4.2f %4.2f %4.2f %4.2f \n", w_col_out[0], w_col_out[1], w_col_out[2], w_col_out[3]);
        sum_cv_value += *cv_value;
        y_col_in = y_col_in + n_row_in;
        w_col_in = w_col_in + n_row_in;
        output_col = output_col + n_row_out;
        w_col_out = w_col_out + n_row_out;
    }
    *cv_value = sum_cv_value;
    return 0;
}

int bin_sparse (double *tt, double *yy, int *ii, int *num_binsp, double *grids,
                double *d, double *data, int *n, double *subject_sum, int *subject_count, double *mintp) {
    int i = 0, j = 0, index = 0, num_bins = *num_binsp, nn = *n;
    double mint = *mintp, dd = *d;
    for (i = 0; i != nn; i++) {
        // may cause problem here
        index = (int) floor((tt[i] - mint) / dd);
        if (index > num_bins) {
            index = num_bins;
        }
        subject_sum[index] += yy[i];
        (subject_count[index])++;
        if (ii[i] != ii[i+1]) {
            for (j = 0; j != num_bins; j++) {
                if (subject_count[j] != 0) {
                    data[num_bins * i + j] = subject_sum[j] / subject_count[j]; 
                    subject_sum[j] = 0; 
                    subject_count[j] = 0; 
                }
            }
        }
    }
    return 0;
}