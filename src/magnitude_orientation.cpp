#include <cmath>
#include <afvec/vectorclass.h>
#include <afvec/vectormath_trig.h>

namespace fasthog {

template <typename Vec_t>
void magnitude_orientation_impl(const double *__restrict gx, const double *__restrict gy, int N, int n_bins,
                                bool signed_hist, double *__restrict magnitude, double *__restrict orientation) {
    const double shift = signed_hist ? 2 * M_PI : M_PI;
    const double scale_factor = n_bins / shift;

    const int n_trunc = Vec_t::size() * (N / Vec_t::size());
    for (int i = 0; i < n_trunc; i += Vec_t::size()) {
        Vec_t GX, GY, MAG, ORIENTATION;
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

template void magnitude_orientation_impl<VEC_T>(const double *__restrict gx, const double *__restrict gy, int N,
                                                int n_bins, bool signed_hist, double *__restrict magnitude,
                                                double *__restrict orientation);

} // namespace fasthog
