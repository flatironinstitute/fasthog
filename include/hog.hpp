#ifndef HOG_HPP
#define HOG_HPP

extern "C" {
void hog(const double *img, int ncols, int nrows, int cell_size_x, int cell_size_y, int block_size_x, int block_size_y,
         int n_bins, double *hist);
void hog_from_gradient(const double *gx, const double *gy, int ncols, int nrows, int cell_size_x, int cell_size_y,
                       int block_size_x, int block_size_y, int n_bins, double *hist);
}
#endif
