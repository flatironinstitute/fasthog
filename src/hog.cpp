#include <hog.hpp>

#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>


#ifdef __x86_64__
#include <cpuid.h>
#include <afvec/vectorclass.h>
#include <afvec/vectormath_trig.h>
#endif

namespace fasthog {
template <typename Vec_t>
void magnitude_orientation_impl(const double *__restrict gx, const double *__restrict gy, int N, int n_bins,
                                bool signed_hist, double *__restrict magnitude, double *__restrict orientation);

constexpr double eps = 1E-5;
constexpr double eps2 = eps * eps;

enum instruction_t : int { Scalar = 0, SSE4 = 1, AVX2 = 2, AVX512 = 3 };

instruction_t get_current_capability() noexcept {
    instruction_t ilvl = Scalar;
#ifdef __x86_64__
    unsigned eax, ebx, ecx, edx, flag = 0;
    int cpuidret = __get_cpuid(1, &eax, &ebx, &ecx, &edx);

    if (!cpuidret)
        return ilvl;

    if (ecx & bit_SSE4_1)
        ilvl = SSE4;

    cpuidret = __get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx);

    if (!cpuidret)
        return ilvl;

    if ((ebx & bit_AVX512F) && (ebx & bit_AVX512VL) && (ebx & bit_AVX512DQ))
        return AVX512;

    if (ebx & bit_AVX2)
        return AVX2;
#endif
    return ilvl;
}

instruction_t get_current_dispatch_target() noexcept {
    const char *dispatch_c_str = std::getenv("FASTHOG_DISPATCH");
    instruction_t user_requested_type = Scalar;

    if (dispatch_c_str) {
        std::string dispatch_str(dispatch_c_str);
        for (auto &c : dispatch_str)
            c = std::tolower(c);

        if (dispatch_str == "scalar") {
            user_requested_type = Scalar;
        } else if (dispatch_str == "sse4") {
            user_requested_type = SSE4;
        } else if (dispatch_str == "avx2") {
            user_requested_type = AVX2;
        } else if (dispatch_str == "avx512") {
            user_requested_type = AVX512;
        }
    }

    instruction_t current_capability = get_current_capability();

    if (dispatch_c_str && user_requested_type > current_capability)
        std::cerr << "WARNING: FASTHOG_DISPATCH environment variable is set to " << dispatch_c_str
                  << ", but the current CPU does not support " << dispatch_c_str << ".\n";

    return dispatch_c_str ? user_requested_type : current_capability;
}
const instruction_t VEC_LEVEL = get_current_dispatch_target();

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

    for (int y = 0; y < ncols; ++y) {
        const int y_cell = y / rows_per_cell;
        const int row_offset = y_cell * n_cells_x * n_bins;
        for (int x = 0; x < nrows; ++x) {
            const int x_cell = x / cols_per_cell;
            const int hist_offset = row_offset + x_cell * n_bins;

            const double angle = orientation[y * ncols + x];
            const double mag = magnitude[y * ncols + x];
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
        }
    }
}

template <>
void magnitude_orientation_impl<double>(const double *__restrict gx, const double *__restrict gy, int N, int n_bins,
                                        bool signed_hist, double *__restrict magnitude,
                                        double *__restrict orientation) {
    const double shift = signed_hist ? 2 * M_PI : M_PI;
    const double scale_factor = n_bins / shift;

    for (int i = 0; i < N; ++i) {
        magnitude[i] = sqrt(gx[i] * gx[i] + gy[i] * gy[i]);
        orientation[i] = atan2(gy[i], gx[i]);
        orientation[i] = scale_factor * (orientation[i] < 0 ? orientation[i] : orientation[i] + shift);
    }
}

void magnitude_orientation(const double *gx, const double *gy, int N, int n_bins, bool signed_hist, double *magnitude,
                           double *orientation) {
#ifdef __x86_64__
    switch (VEC_LEVEL) {
    case Scalar:
        magnitude_orientation_impl<double>(gx, gy, N, n_bins, signed_hist, magnitude, orientation);
        return;
    case SSE4:
        magnitude_orientation_impl<Vec2d>(gx, gy, N, n_bins, signed_hist, magnitude, orientation);
        return;
    case AVX2:
        magnitude_orientation_impl<Vec4d>(gx, gy, N, n_bins, signed_hist, magnitude, orientation);
        return;
    case AVX512:
        magnitude_orientation_impl<Vec8d>(gx, gy, N, n_bins, signed_hist, magnitude, orientation);
        return;
    }
#else
    magnitude_orientation_impl<double>(gx, gy, N, n_bins, signed_hist, magnitude, orientation);
#endif
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
} // namespace fasthog

extern "C" {
void fasthog_hog(const double *img, int ncols, int nrows, int cell_size_x, int cell_size_y, int block_size_x,
                 int block_size_y, int n_bins, bool signed_hist, NORM_TYPE norm_type, double *hist) {
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
    using namespace fasthog;

    gradient(img, nrows, ncols, gx, gy);
    magnitude_orientation(gx, gy, N_pixels, n_bins, signed_hist, magnitude, orientation);
    build_histogram(magnitude, orientation, nrows, ncols, cell_size_y, cell_size_x, n_bins, unblocked_hist);

    normalize_histogram(unblocked_hist, n_cells_x, n_cells_y, block_size_x, block_size_y, n_bins, norm_type, hist);

    delete[] mempool;
}

void fasthog_hog_from_gradient(const double *gx, const double *gy, int ncols, int nrows, int cell_size_x,
                               int cell_size_y, int block_size_x, int block_size_y, int n_bins, bool signed_hist,
                               NORM_TYPE norm_type, double *hist) {
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
    using namespace fasthog;

    magnitude_orientation(gx, gy, N_pixels, n_bins, signed_hist, magnitude, orientation);
    build_histogram(magnitude, orientation, nrows, ncols, cell_size_y, cell_size_x, n_bins, unblocked_hist);
    normalize_histogram(unblocked_hist, n_cells_x, n_cells_y, block_size_x, block_size_y, n_bins, norm_type, hist);

    delete[] mempool;
}
}
