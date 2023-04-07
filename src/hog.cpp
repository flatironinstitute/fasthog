#include <cstring>
#include <hog.hpp>

#include <cmath>
#include <memory>
#include <cstdio>

#include <afvec/vectorclass.h>
#include <afvec/vectormath_trig.h>

static const double eps = 1E-7;

extern "C" {
void build_histogram(const double *magnitude, const double *orientation, int nrows, int ncols, int rows_per_cell,
                     int cols_per_cell, int n_bins, double *res) {
    const int n_cells_y = nrows / rows_per_cell;
    const int n_cells_x = ncols / cols_per_cell;

    const double half_bin = 180.0 / n_bins;
    const double inv_bin_size = n_bins / 360.0;
    double bin_centers[n_bins + 2];
    bin_centers[0] = -half_bin;
    for (int i = 1; i < n_bins + 2; ++i)
        bin_centers[i] = bin_centers[i - 1] + 2 * half_bin;

    for (int y_cell = 0; y_cell < n_cells_y; ++y_cell) {
        const int y_start = y_cell * rows_per_cell;
        const int y_end = y_start + rows_per_cell;
        for (int x_cell = 0; x_cell < n_cells_x; ++x_cell) {
            const int x_start = x_cell * cols_per_cell;
            const int x_end = x_start + cols_per_cell;
            const int hist_offset = y_cell * n_cells_x * n_bins + x_cell * n_bins;

            double norm_factor = 0.0;
            for (int y = y_start; y < y_end; ++y) {
                for (int x = x_start; x < x_end; ++x) {
                    double angle = orientation[y*ncols + x];
                    double mag = magnitude[y*ncols + x];
                    int low_bin = inv_bin_size * (angle - half_bin);
                    int high_bin = low_bin + 1;
                    double low_center = bin_centers[low_bin + 1];
                    double high_center = bin_centers[high_bin + 1];
                    if (angle - half_bin < 0) {
                        low_bin = n_bins - 1;
                        low_center = bin_centers[0];
                    }
                    if (high_bin > n_bins - 1)
                        high_bin = 0;

                    res[hist_offset + low_bin] = mag * inv_bin_size * (angle - low_center);
                    res[hist_offset + high_bin] = mag * inv_bin_size * (high_center - angle);
                    norm_factor += mag;
                }
            }

            norm_factor = 1.0 / sqrt(norm_factor * norm_factor + eps);
            double norm_factor2 = 0.0;
            for (int i_bin = 0; i_bin < n_bins; ++i_bin)  {
                res[hist_offset + i_bin] *= norm_factor;
                res[hist_offset + i_bin] = std::min(0.2, res[hist_offset + i_bin]);
                norm_factor2 += res[hist_offset + i_bin];
            }

            norm_factor2 = 1.0 / sqrt(norm_factor2 * norm_factor2 + eps);
            for (int i_bin = 0; i_bin < n_bins; ++i_bin)
                res[hist_offset + i_bin] *= norm_factor2;
        }
    }
}

void magnitude_orientation(const double *gx, const double *gy, int N, double *magnitude, double *orientation) {
    const double scale = 180.0 / M_PI;
    for (int i = 0; i < N; i += 4) {
        Vec4d GX, GY, MAG, ORIENTATION;
        GX.load(gx + i);
        GY.load(gy + i);
        MAG = sqrt(GX * GX + GY * GY);
        ORIENTATION = atan2(GY, GX) * scale;
        ORIENTATION = if_add(ORIENTATION < 0, ORIENTATION, 360.0);
        MAG.store(magnitude + i);
        ORIENTATION.store(orientation + i);
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

void hog(const double *img, int ncols, int nrows, int cell_size_x, int cell_size_y, int n_bins, double *hist) {
    const int N = nrows * ncols;
    const int n_cells_x = ncols / cell_size_x;
    const int n_cells_y = nrows / cell_size_x;
    double gx[N], gy[N], magnitude[N], orientation[N];
    gradient(img, nrows, ncols, gx, gy);
    magnitude_orientation(gx, gy, N, magnitude, orientation);
    memset(hist, 0, n_cells_x * n_cells_y * n_bins * sizeof(double));
    build_histogram(magnitude, orientation, nrows, ncols, cell_size_y, cell_size_x, n_bins, hist);
}
}
