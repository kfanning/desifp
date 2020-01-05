import os
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder

filename = 'fvc.20200104204051.fits'
path = os.path.join('/data/images/fvc', filename)
print(f'Detecting sources in: {path}')
with fits.open(path) as hdul:
    image = hdul[0].data
mean, median, std = sigma_clipped_stats(image, sigma=3)
daofind = DAOStarFinder(fwhm=3, threshold=100*std)
sources = daofind(image - median)
for col in sources.colnames:
    sources[col].info.format = '%.8g'
print(sources)
