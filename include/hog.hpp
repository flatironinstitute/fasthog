#ifndef HOG_HPP
#define HOG_HPP

extern "C" {
enum NORM_TYPE {
    NONE = 0,
    L1 = 1,
    L1_SQRT = 2,
    L2 = 3,
    L2_HYS = 4,
};
void hog(const double *img, int ncols, int nrows, int cell_size_x, int cell_size_y, int block_size_x, int block_size_y,
         int n_bins, bool signed_hist, NORM_TYPE norm_type, double *hist);
void hog_from_gradient(const double *gx, const double *gy, int ncols, int nrows, int cell_size_x, int cell_size_y,
                       int block_size_x, int block_size_y, int n_bins, bool signed_hist, NORM_TYPE norm_type,
                       double *hist);
}
#endif
