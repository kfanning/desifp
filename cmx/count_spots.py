import os
import glob
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder, CircularAperture, aperture_photometry
import matplotlib.pyplot as plt

# filename = 'fvc.20200104204051.fits'
# find the latest fvc image
paths = glob.glob('/data/images/fvc/*.fits')
path_latest = max(paths, key=os.path.getctime)
print(f'Detecting sources in: {path_latest}')
with fits.open(path_latest) as hdul:
    image = hdul[0].data
mean, median, std = sigma_clipped_stats(image, sigma=3)
daofind = DAOStarFinder(fwhm=3, threshold=100*std)
sources = daofind(image - median)
for col in sources.colnames:
    sources[col].info.format = '%.8g'
sources.sort('peak')
print(sources)
xc = sources['xcentroid'].data
yc = sources['ycentroid'].data
aperture = CircularAperture(np.array([xc, yc]).T, r=4)
phot = aperture_photometry(image-median, aperture)
print(f'{len(sources)} sources detected in: {path_latest}')
spectcon_num = input('Enter the 1-digit enabled SPECTCON#: ')
duty = input('Enter the duty cycle used (e.g. 0.10): ')
if spectcon_num and duty:
    n_bins = 30
    fig, ax = plt.subplots()
    ax.hist(sources['peak'], bins=n_bins, color='C0', histtype='step',
            label='Peak ADU')
    ax.set_xlabel('Peak ADU')
    ax.set_ylabel('Count of sources')
    ax.legend(loc=1)
    # other plot
    ax = ax.twiny()
    ax.hist(phot['aperture_sum'], bins=n_bins, color='C1', histtype='step',
            label='Aperture Sum')
    ax.set_xlabel('Aperture Sum')

    ax.legend(loc=2)
    ax.set_title(f'FVC image illumination distribution:'
                 f'\nSPECTCON{spectcon_num}, duty_cycle = {duty}'
                 f'\n{os.path.basename(path_latest)}')
    fig.savefig(
        f'fvc_illumination_dist-SPECTCON{spectcon_num}-duty_{duty}.pdf',
        bbox_inches='tight')
