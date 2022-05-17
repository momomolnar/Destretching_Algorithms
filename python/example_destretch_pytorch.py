#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 21:23:57 2020

@author: molnarad
"""
from astropy.io import fits

from destretch_pytorch import *

from flicker import flicker


Hdu = fits.open("test/test.fits")

test_data = Hdu[0].data

scene = np.array(test_data[4, :, :]).astype(np.float32)
reference = np.array(test_data[10, :, :]).astype(np.float32)

kernel_sizes = [100, 64, 32]


test_destretch(scene, reference, kernel_sizes,
               plot=True, device="GPU")
