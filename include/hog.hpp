#ifndef HOG_HPP
#define HOG_HPP

extern "C" {
void hog(const double *img, int ncols, int nrows, int cell_size_x, int cell_size_y, int n_bins, double *hist);
}
#endif
