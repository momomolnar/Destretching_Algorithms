#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 21:23:57 2020

@author: molnarad
"""

from destretch import * 

from astropy.io import fits

hdu = fits.open("test/test.fits")

test_data = hdu[0].data

scene        = test_data[4, :, :]
reference    = test_data[10, :, :] 
kernel_sizes = [100, 64, 32]

test_destretch(scene, reference, kernel_sizes)
