# FastHOG
Reasonably fast implementation of Histogram of Oriented Gradients. Designed to be roughly equivalent to
`skimage.feature.hog`. While slightly less flexible (it only supports `float64` single-channel images),
it is _significantly_ faster, so ideal in workloads where `HOG` calculation is a bottleneck, as was in
the https://github.com/flatironinstitute/ManifoldEM project that inspired this repository.


## Installation

```bash
pip install fasthog
```

## Usage
Example taken from https://scikit-image.org/docs/stable/auto_examples/features_detection/plot_hog.html

```python3
from fasthog import hog
from skimage import data
from skimage.color import rgb2gray

image = rgb2gray(data.astronaut())

fd, hog_image = hog(
    image,
    n_bins=8,
    pixels_per_cell=(16, 16),
    cells_per_block=(1, 1),
)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), sharex=True, sharey=True)

ax1.axis('off')
ax1.imshow(image, cmap=plt.cm.gray)
ax1.set_title('Input image')

# Rescale histogram for better display
hog_image_rescaled = exposure.rescale_intensity(hog_image, in_range=(0, 10))

ax2.axis('off')
ax2.imshow(hog_image_rescaled, cmap=plt.cm.gray)
ax2.set_title('Histogram of Oriented Gradients')
plt.show()


```
