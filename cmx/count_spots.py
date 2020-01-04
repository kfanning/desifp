import os
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder

filename = 'fvc.00037412_0000.fits'
path = os.path.join('/data/images/fvc', filename)
with fits.open(path) as hdul:
    image = hdul[0].data
mean, median, std = sigma_clipped_stats(image, sigma=3)
daofind = DAOStarFinder(fwhm=3, threshold=200*std)
sources = daofind(image - median)
for col in sources.colnames:
    sources[col].info.format = '%.8g'
print(sources)
