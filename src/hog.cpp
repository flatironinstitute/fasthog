#include <hog.hpp>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>

#include <afvec/vectorclass.h>
#include <afvec/vectormath_trig.h>

constexpr double eps = 1E-5;
constexpr double eps2 = eps * eps;

void normalize_histogram(const double *unblocked, int n_cells_x, int n_cells_y, int block_size_x, int block_size_y,
                         int n_bins, NORM_TYPE norm_type, double *__restrict__ hist) {
    const int n_blocks_x = (n_cells_x - block_size_x + 1);
    const int n_blocks_y = (n_cells_y - block_size_y + 1);
    memset(hist, 0, n_blocks_x * n_blocks_y * n_bins * sizeof(double));

    for (int y_block = 0; y_block < n_blocks_y; ++y_block) {
        for (int x_block = 0; x_block < n_blocks_x; ++x_block) {
            double *block = hist + y_block * n_blocks_x * n_bins + x_block * n_bins;
            for (int y_cell = y_block; y_cell < y_block + block_size_y; ++y_cell) {
                for (int x_cell = x_block; x_cell < x_block + block_size_x; ++x_cell) {
                    const double *cell = unblocked + (y_cell * n_cells_x + x_cell) * n_bins;
                    for (int i_bin = 0; i_bin < n_bins; ++i_bin)
                        block[i_bin] += cell[i_bin];
                }
            }
        }
    }

    switch (norm_type) {
    case NONE: {
        break;
    }
    case L1: {
        for (int i_block = 0; i_block < n_blocks_x * n_blocks_y; ++i_block) {
            double *block = hist + i_block * n_bins;
            double norm_factor = eps;
            for (int i_bin = 0; i_bin < n_bins; ++i_bin)
                norm_factor += block[i_bin];

            norm_factor = 1.0 / norm_factor;
            for (int i_bin = 0; i_bin < n_bins; ++i_bin)
                block[i_bin] *= norm_factor;
        }
        break;
    }
    case L1_SQRT: {
        for (int i_block = 0; i_block < n_blocks_x * n_blocks_y; ++i_block) {
            double *block = hist + i_block * n_bins;
            double norm_factor = eps;
            for (int i_bin = 0; i_bin < n_bins; ++i_bin)
                norm_factor += block[i_bin];

            norm_factor = 1.0 / norm_factor;
            for (int i_bin = 0; i_bin < n_bins; ++i_bin)
                block[i_bin] = sqrt(norm_factor * block[i_bin]);
        }
        break;
    }
    case L2: {
        for (int i_block = 0; i_block < n_blocks_x * n_blocks_y; ++i_block) {
            double *block = hist + i_block * n_bins;
            double norm_factor = eps2;
            for (int i_bin = 0; i_bin < n_bins; ++i_bin)
                norm_factor += block[i_bin] * block[i_bin];

            norm_factor = 1.0 / sqrt(norm_factor);
            for (int i_bin = 0; i_bin < n_bins; ++i_bin)
                block[i_bin] *= norm_factor;
        }
        break;
    }
    case L2_HYS: {
        for (int i_block = 0; i_block < n_blocks_x * n_blocks_y; ++i_block) {
            double *block = hist + i_block * n_bins;
            double norm_factor = eps2;
            for (int i_bin = 0; i_bin < n_bins; ++i_bin)
                norm_factor += block[i_bin] * block[i_bin];

            norm_factor = 1.0 / sqrt(norm_factor);

            double norm_factor2 = eps2;
            for (int i_bin = 0; i_bin < n_bins; ++i_bin) {
                block[i_bin] *= norm_factor;
                block[i_bin] = std::min(0.2, block[i_bin]);
                norm_factor2 += block[i_bin] * block[i_bin];
            }

            norm_factor2 = 1.0 / sqrt(norm_factor2);
            for (int i_bin = 0; i_bin < n_bins; ++i_bin)
                block[i_bin] *= norm_factor2;
        }
        break;
    }
    }
}

void build_histogram(const double *magnitude, const double *orientation, int nrows, int ncols, int rows_per_cell,
                     int cols_per_cell, int n_bins, double *hist) {
    const int n_cells_y = nrows / rows_per_cell;
    const int n_cells_x = ncols / cols_per_cell;
    memset(hist, 0, n_cells_x * n_cells_y * n_bins * sizeof(double));
    constexpr bool interp = false;

    for (int y = 0; y < ncols; ++y) {
        const int y_cell = y / rows_per_cell;
        const int row_offset = y_cell * n_cells_x * n_bins;
        for (int x = 0; x < nrows; ++x) {
            const int x_cell = x / cols_per_cell;
            const int hist_offset = row_offset + x_cell * n_bins;

            const double angle = orientation[y * ncols + x];
            const double mag = magnitude[y * ncols + x];
            if (interp) {
                int high_bin = angle + 0.5;
                int low_bin = high_bin - 1;

                const double low_vote = mag * (high_bin + 0.5 - angle);
                const double high_vote = mag - low_vote;
                if (high_bin < 1)
                    low_bin = n_bins - 1;
                if (high_bin >= n_bins)
                    high_bin = 0;

                hist[hist_offset + low_bin] += low_vote;
                hist[hist_offset + high_bin] += high_vote;
            } else {
                int bin = angle;
                hist[hist_offset + bin] += mag;
            }
        }
    }
}

void magnitude_orientation(const double *gx, const double *gy, int N, int n_bins, bool signed_hist, double *magnitude,
                           double *orientation) {
    double shift = signed_hist ? 2 * M_PI : M_PI;
    double scale_factor = n_bins / shift;

    int n_trunc = 4 * (N / 4);
    for (int i = 0; i < n_trunc; i += 4) {
        Vec4d GX, GY, MAG, ORIENTATION;
        GX.load(gx + i);
        GY.load(gy + i);
        MAG = sqrt(GX * GX + GY * GY);
        ORIENTATION = atan2(GY, GX);
        ORIENTATION = scale_factor * if_add(ORIENTATION < 0, ORIENTATION, shift);
        MAG.store(magnitude + i);
        ORIENTATION.store(orientation + i);
    }

    for (int i = n_trunc; i < N; ++i) {
        magnitude[i] = sqrt(gx[i] * gx[i] + gy[i] * gy[i]);
        orientation[i] = atan2(gy[i], gx[i]);
        orientation[i] = scale_factor * (orientation[i] < 0 ? orientation[i] : orientation[i] + shift);
    }
}

void gradient(const double *img, int nrows, int ncols, double *gx, double *gy) {
    for (int y = 0; y < nrows; ++y) {
        gx[y * ncols] = -img[y * ncols] + img[y * ncols + 1];
        const int yoff = y * ncols;
        for (int x = 1; x < ncols - 1; ++x)
            gx[yoff + x] = -img[yoff + x - 1] + img[yoff + x + 1];
        gx[(y + 1) * ncols - 1] = -img[(y + 1) * ncols - 2] + img[(y + 1) * ncols - 1];
    }

    for (int x = 0; x < ncols; ++x) {
        gy[x] = img[x] - img[ncols + x];
        for (int y = 1; y < nrows - 1; ++y)
            gy[y * ncols + x] = img[(y - 1) * ncols + x] - img[(y + 1) * ncols + x];
        gy[(nrows - 1) * ncols + x] = img[(nrows - 2) * ncols + x] - img[(nrows - 1) * ncols + x];
    }
}

extern "C" {
void hog(const double *img, int ncols, int nrows, int cell_size_x, int cell_size_y, int block_size_x, int block_size_y,
         int n_bins, bool signed_hist, NORM_TYPE norm_type, double *hist) {
    const int N_pixels = nrows * ncols;
    const int n_cells_x = ncols / cell_size_x;
    const int n_cells_y = nrows / cell_size_y;
    const int N_cells = n_cells_x * n_cells_y;
    const int n_blocks_x = (n_cells_x - block_size_x) + 1;
    const int n_blocks_y = (n_cells_y - block_size_y) + 1;
    double *mempool = new double[4 * N_pixels + N_cells * n_bins];
    double *gx = mempool + 0 * N_pixels;
    double *gy = mempool + 1 * N_pixels;
    double *magnitude = mempool + 2 * N_pixels;
    double *orientation = mempool + 3 * N_pixels;
    double *unblocked_hist = mempool + 4 * N_pixels;

    gradient(img, nrows, ncols, gx, gy);
    magnitude_orientation(gx, gy, N_pixels, n_bins, signed_hist, magnitude, orientation);
    build_histogram(magnitude, orientation, nrows, ncols, cell_size_y, cell_size_x, n_bins, unblocked_hist);

    normalize_histogram(unblocked_hist, n_cells_x, n_cells_y, block_size_x, block_size_y, n_bins, norm_type, hist);

    delete[] mempool;
}

void hog_from_gradient(const double *gx, const double *gy, int ncols, int nrows, int cell_size_x, int cell_size_y,
                       int block_size_x, int block_size_y, int n_bins, bool signed_hist, NORM_TYPE norm_type,
                       double *hist) {
    const int N_pixels = nrows * ncols;
    const int n_cells_x = ncols / cell_size_x;
    const int n_cells_y = nrows / cell_size_y;
    const int N_cells = n_cells_x * n_cells_y;
    const int n_blocks_x = (n_cells_x - block_size_x) + 1;
    const int n_blocks_y = (n_cells_y - block_size_y) + 1;
    double *mempool = new double[2 * N_pixels + N_cells * n_bins];
    double *magnitude = mempool + 0 * N_pixels;
    double *orientation = mempool + 1 * N_pixels;
    double *unblocked_hist = mempool + 2 * N_pixels;

    magnitude_orientation(gx, gy, N_pixels, n_bins, signed_hist, magnitude, orientation);
    build_histogram(magnitude, orientation, nrows, ncols, cell_size_y, cell_size_x, n_bins, unblocked_hist);
    normalize_histogram(unblocked_hist, n_cells_x, n_cells_y, block_size_x, block_size_y, n_bins, norm_type, hist);

    delete[] mempool;
}
}
