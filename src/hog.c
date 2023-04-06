#include <hog.h>

#include <cmath>

void magnitude_orientation(const double *gx, const double *gy, int N, double *magnitude, double *orientation) {
    const double scale = 180.0 / M_PI;
    for (int i = 0; i < N; ++i) {
        magnitude[i] = sqrt(gx[i] * gx[i] + gy[i] * gy[i]);
        orientation[i] = atan2(gy[i], gx[i]) * scale;
    }
}

void gradient(const double *img, int ny, int nx, double *gx, double *gy) {
    for (int y = 0; y < ny; ++y) {
        gx[y * nx] = -img[y * nx] + img[y * nx + 1];
        const int yoff = y * nx;
        for (int x = 1; x < nx - 1; ++x)
            gx[yoff + x] = -img[yoff + x - 1] + img[yoff + x + 1];
        gx[(y + 1) * nx - 1] = -img[(y + 1) * nx - 2] + img[(y + 1) * nx - 1];
    }

    for (int x = 0; x < nx; ++x) {
        gy[x] = img[x] - img[nx + x];
        for (int y = 1; y < ny - 1; ++y)
            gy[y * nx + x] = img[(y - 1) * nx + x] - img[(y + 1) * nx + x];
        gy[(ny - 1) * nx + x] = img[(ny - 2) * nx + x] - img[(ny - 1) * nx + x];
    }
}

void hog(double *img, int ny, int nx) {
    const int N = ny * nx;
    double gx[N], gy[N], magnitude[N], orientation[N];
    gradient(img, ny, nx, gx, gy);
    magnitude_orientation(gx, gy, N, magnitude, orientation);
}
