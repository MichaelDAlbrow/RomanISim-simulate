"""
Convert an asdf image file to fits format.

Basic code from Andrea Bellini via Matthew Penny.

Michael Albrow
"""

import sys
import numpy as np
import astropy
from astropy.io import fits
import roman_datamodels as rdm


def asdf_to_fits(asdf_file: str, fits_file: (str, type(None)) = None) -> None:

    imagehandler = rdm.open(asdf_file)
    image = np.array(imagehandler.data)
    err = np.array(imagehandler.err)
    dq = np.array(imagehandler.dq)

    # header info
    print(imagehandler.meta)

    # saving the wcs info
    wcs_info = imagehandler.meta.wcs.to_fits(bounding_box=((0, 4095), (0, 4095)))

    # Converting the output into a ds9 readable format for quick visual inspection
    data_HDU = fits.PrimaryHDU(image, header=wcs_info[0])
    err_HDU = fits.ImageHDU(data=err, name="ERR")
    dq_HDU = fits.ImageHDU(data=dq, name="DQ")
    HDUl = fits.HDUList([data_HDU, err_HDU, dq_HDU])

    if fits_file is None:
        fits_file = asdf_file.split('.')[0] + '.fits'

    HDUl.writeto(fits_file, overwrite=True)


if __name__ == '__main__':

    asdf_file = sys.argv[1]

    fits_file = None
    try:
        fits_file = sys.argv[2]
    except IndexError:
        pass

    asdf_to_fits(asdf_file, fits_file)

