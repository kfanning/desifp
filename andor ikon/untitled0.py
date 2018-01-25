# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 02:37:33 2016

@author: givoltage
"""

#%% import libraries
import sys
print(sys.prefix)

import os
import numpy as np
import matplotlib.pyplot as plt
# import cv2
from astropy.io import fits
from skimage import measure

#%% read in file
dir = r"C:\Users\givoltage\Google Drive\DESI\protoDESI\images"
os.chdir(dir)
hdulist = fits.open("5W_10fib_lens_0s_001-001.fits")
header = hdulist[0].header
img = hdulist[0].data
#data = hdulist[0].data[::-1, :]

# creating binary mask

img_binary = img > intensity_threshold
img_renorm = img/np.amax(img)
plt.imshow(img_binary, origin='lower')
img_label = measure.label(img_binary, connectivity=img.ndim)
plt.imshow(img_label, origin='lower')
# getting region properties
# for list of available properties, check documentation below
# http://scikit-image.org/docs/dev/api/skimage.measure.html?highlight=ellipse#regionprops

# calculate contours and plot boundaries
contours = measure.find_contours(img, intensity_threshold)
[fit, ax] = plt.subplots(figsize=(3.354, 2.529),dpi=100)
ax.imshow(img,interpolation='nearest', cmap=plt.cm.gray,
										origin='lower')
# consistency check
if len(contours) == len(imgProps):
	for i, contour in enumerate(contours):
		ax.plot(contour[:, 1], contour[:, 0], linewidth=0.1, color = 'b')
	for i, imgProp in enumerate(imgProps):
		ax.plot(imgProp.weighted_centroid[1], imgProp.weighted_centroid[0],'rx', markersize = 0.1)
	ax.axis('image')
	#ax.set_xticks([])
	#ax.set_yticks([])
	mydpi = 1200
	plt.savefig('5W_10fib_lens_0s_001-001.png', dpi=mydpi)
	#plt.savefig('1_001.png', dpi=mydpi, bbox_inches='tight', pad_inches=0)
	# plt.savefig('1_001.svg', dpi=mydpi)
	# plt.savefig('1_001.pdf', dpi=mydpi)
else:
	print("dimensions don't match between contours and properties")





### test section
hdu = fits.open('dataset.fits')[0]

### test section
from skimage import data, util
from skimage.measure import label
coin = data.coins()
coinmask = util.img_as_ubyte(coin) > 110
label_coin = label(coinmask, connectivity=img.ndim)
plt.imshow(coinmask, origin='lower')
props = measure.regionprops(label_coin)
# centroid of first labeled object
props[0].centroid










# Construct some test data
[x, y] = np.ogrid[-np.pi:np.pi:100j, -np.pi:np.pi:100j]
r = np.sin(np.exp((np.sin(x)**3 + np.cos(y)**2)))

# Find contours at a constant value of 0.8
contours = measure.find_contours(r, 0.8)

# Display the image and plot all contours found
[fig, ax] = plt.subplots(figsize=(10, 10))
ax.imshow(r, interpolation='nearest', cmap=plt.cm.gray)

for n, contour in enumerate(contours):
    ax.plot(contour[:, 1], contour[:, 0], linewidth=2m, colors='k')

ax.axis('image')
ax.set_xticks([])
ax.set_yticks([])
plt.show()





---------------------------------------------------
# creating binary mask
intensity_threshold = 2100
img_binary = img > intensity_threshold
img_renorm = img/np.amax(img)
plt.imshow(img_binary, origin='lower')
img_label = measure.label(img_binary, connectivity=img.ndim)
plt.imshow(img_label, origin='lower')
# getting region properties
# for list of available properties, check documentation below
# http://scikit-image.org/docs/dev/api/skimage.measure.html?highlight=ellipse#regionprops

# calculate contours and plot boundaries
contours = measure.find_contours(img, intensity_threshold)
[fit, ax] = plt.subplots(figsize=(3.354, 2.529),dpi=100)
ax.imshow(img,interpolation='nearest', cmap=plt.cm.gray,
										origin='lower')
# consistency check
if len(contours) == len(imgProps):
	for i, contour in enumerate(contours):
		ax.plot(contour[:, 1], contour[:, 0], linewidth=0.1, color = 'b')
	for i, imgProp in enumerate(imgProps):
		ax.plot(imgProp.weighted_centroid[1], imgProp.weighted_centroid[0],'rx', markersize = 0.1)
	ax.axis('image')
	#ax.set_xticks([])
	#ax.set_yticks([])
	mydpi = 1200
	plt.savefig('5W_10fib_lens_0s_001-001.png', dpi=mydpi)
	#plt.savefig('1_001.png', dpi=mydpi, bbox_inches='tight', pad_inches=0)
	# plt.savefig('1_001.svg', dpi=mydpi)
	# plt.savefig('1_001.pdf', dpi=mydpi)
else:
	print("dimensions don't match between contours and properties")





### test segmentation
import numpy as np
from photutils import datasets
hdu = datasets.load_star_image()
data = hdu.data[0:400, 0:400]
image = hdu.data.astype(float)
image -= np.median(image)

from photutils import daofind
from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats
bkg_sigma = mad_std(image)
mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5)
print(mean, median, std)
sources = daofind(image, fwhm=4.0, threshold=3.0*bkg_sigma)
print(sources)

from photutils import CircularAperture
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import matplotlib.pylab as plt
positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(positions, r=4.)
norm = ImageNormalize(stretch=SqrtStretch())
plt.imshow(data, cmap='Greys', origin='lower', norm=norm)
apertures.plot(color='blue', lw=1.5, alpha=0.5)

#
from photutils.datasets import make_100gaussians_image
data = make_100gaussians_image()
from photutils import detect_threshold
threshold = detect_threshold(data, snr=3.)

