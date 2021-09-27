#!/usr/bin/env python3

from destretch import *
import time
import matplotlib.pyplot as plt
from scipy.ndimage import rotate, shift, gaussian_filter

def radial_distances(array_size, central_point=[0,0]):

    if len(array_size) > 2:
        sx = array_size[1]
        sy = array_size[2]
    elif len(array_size) > 1:
        sx = array_size[0]
        sy = array_size[1]
    else:
        sx = array_size[0]
        sy = array_size[0]
    if not np.any(central_point):
        central_point = [sx/2., sy/2.]
    xdist = np.tile(np.arange(sx) - central_point[0], (sy, 1))
    ydist = np.tile(np.arange(sy).reshape((-1, 1)) - central_point[1], (1, sx))
    radist = np.sqrt(xdist**2 + ydist**2)
    
    return radist


def rotate_around_point(img, angle, pivot, interp=3):

    padX = [img.shape[1] - pivot[0], pivot[0]]
    padY = [img.shape[0] - pivot[1], pivot[1]]
    imgP = np.pad(img, [padY, padX], 'constant')
    imgR = rotate(imgP, angle, reshape=False, order=interp)
    return imgR[padY[0] : -padY[1], padX[0] : -padX[1]]


def calc_synthetic_destretch(kernel_size,
                             apod_percent=0.08,
                             input_image=0,
                             output_image=0,
                             regenerate_image=0,
                             dot_radius=5.1,
                             dot_step=12,
                             apply_rotation=0,
                             apply_shift=[0, 0],
                             add_noise=0,
                             noise_level=0,
                             debug_level=0,
                             plot_level=0,
                             flip_image=0,
                             gaussian_blur=0.0):

    rot_angle = apply_rotation
    cent = 500
    spacrat = 0.5
    sft_x = apply_shift[0]
    sft_y = apply_shift[1]
    buffer = 25
    
    if (not input_image) or regenerate_image:
        print('Generating random images')
        # Generate image of random dots
        numx = 1000
        numy = 1000
        nstep = int(numx / dot_step)
        dot_loc_randomness = 2.0
        dots_image_ref = np.zeros((numx + buffer*2, numy + buffer*2), order='F')
        dots_image_sft = np.zeros((numx + buffer*2, numy + buffer*2), order='F')
        xcoord = np.arange(nstep) * numx/nstep
        ycoord = np.arange(nstep) * numy/nstep
        xcoord += (numx - xcoord.max()) / 2.
        ycoord += (numy - ycoord.max()) / 2.
        rng = np.random.default_rng()
        coord_rand = rng.normal(size=(2, nstep, nstep)) * dot_step/dot_loc_randomness
        coord_rand = np.array(coord_rand, order='F')
        xcoord_num = xcoord.size
        ycoord_num = ycoord.size
        dot_coord_ref = np.zeros((2, xcoord_num, ycoord_num), order='F')
        
        for x in range(xcoord_num):
            for y in range(ycoord_num):
                xcoord_jitter = np.clip(xcoord[x] + coord_rand[0, x, y], 0, numx-1)
                ycoord_jitter = np.clip(ycoord[y] + coord_rand[1, x, y], 0, numy-1)
                dot_coord_ref[:, x, y] = [xcoord_jitter, ycoord_jitter]
        
        if apply_rotation:
            print('Rotating random images')
            beta = np.arctan2(dot_coord_ref[0] - cent, dot_coord_ref[1] - cent)
            alpha = rot_angle * np.pi/180.
            distance = np.sqrt(np.sum((dot_coord_ref - cent)**2, 0))
            shifts_calc_x = np.sin(beta + alpha) * distance - np.sin(beta) * distance
            shifts_calc_y = np.cos(beta + alpha) * distance - np.cos(beta) * distance
        
        for x in range(xcoord_num):
            for y in range(ycoord_num):
                if apply_rotation:
                    offsets = [shifts_calc_x[y, x], shifts_calc_y[y, x]]
                else:
                    offsets = apply_shift
                subfield_pix = np.around(dot_coord_ref[:, x, y]).astype(int) - buffer
                dot_subfield = np.clip(dot_radius - radial_distances([buffer*2], dot_coord_ref[:, x, y] - subfield_pix), 0, None)
                dots_image_ref[subfield_pix[0] + buffer:subfield_pix[0] + buffer*3, subfield_pix[1] + buffer:subfield_pix[1] + buffer*3] += dot_subfield
                dot_subfield = np.clip(dot_radius - radial_distances([buffer*2], dot_coord_ref[:, x, y] + offsets - subfield_pix), 0, None)
                dots_image_sft[subfield_pix[0] + buffer:subfield_pix[0] + buffer*3, subfield_pix[1] + buffer:subfield_pix[1] + buffer*3] += dot_subfield
        dots_image_ref = dots_image_ref[buffer:buffer+numx, buffer:buffer+numy]
        dots_image_sft = dots_image_sft[buffer:buffer+numx, buffer:buffer+numy]
        print('Images successfully created')
    else:
        print('Using provided images')
        dots_image_ref = input_image
        scene_size = dots_image_ref.shape
        numx = scene_size[0]
        numy = scene_size[1]
        
        if apply_rotation:
            dots_image_sft = rotate_around_point(dots_image_ref, rot_angle, [cent, cent], 3)
        else:
            dots_image_sft = shift(dots_image_ref, apply_shift)
    
#    corr = signal.correlate2d(dots_image_sft, dots_image_ref, boundary='symm')
    corr = np.fft.ifft2(np.fft.fft2(dots_image_sft).conj()*np.fft.fft2(dots_image_ref))
    calc_shift = np.array(np.unravel_index(np.argmax(corr), corr.shape))
    print('calc_shift = ' + repr(calc_shift))
    if gaussian_blur:
        print('Degrading with gaussian blur')
        dots_image_ref = gaussian_filter(dots_image_ref, gaussian_blur)
        dots_image_sft = gaussian_filter(dots_image_sft, gaussian_blur)
    else:
        print('No blurring')
    if plot_level:
        print('Plotting reference image')
        plt.ion()
        plt.imshow(dots_image_ref)
    else:
        print('No plotting')
    
    if add_noise:
        print('Degrading image with shot noise')
        rng2 = np.random.default_rng()
        rng3 = np.random.default_rng()
        data_ref_norm = dots_image_ref / np.mean(dots_image_ref)
        data_sft_norm = dots_image_sft / np.mean(dots_image_sft)
        dots_image_ref = data_ref_norm * (1 + rng2.normal(size=(numx, numy))*noise_level)
        dots_image_sft = data_sft_norm * (1 + rng3.normal(size=(numx, numy))*noise_level)
    else:
        print('No shot noise')
    
    if flip_image:
        print('Transposing images')
        dots_image_ref = dots_image_ref.T
        dots_image_sft = dots_image_sft.T
    output_image = dots_image_ref
    
    print('Now destretching shifted image, please wait...')
    data_sft_noisy_reg, disp, rdisp, d_info = reg(dots_image_sft, dots_image_ref, kernel_size)
    print('Destretching complete!')
    shifts_destr_x = (disp - rdisp)[0]
    shifts_destr_y = (disp - rdisp)[1]
    
    if apply_rotation:
        beta = np.arctan2(rdisp[0] - cent, rdisp[1] - cent)
        alpha = rot_angle * np.pi/180.
        distance = np.sqrt(np.sum((rdisp - cent)**2, 0))
        shifts_calc_x = np.sin(beta + alpha) * distance - np.sin(beta) * distance
        shifts_calc_y = np.cos(beta + alpha) * distance - np.cos(beta) * distance
        vector_ratio = [np.median(shifts_destr_x/shifts_calc_x), np.median(shifts_destr_y/shifts_calc_y)]
    else:
        shifts_calc_x = rdisp[0]*0 + sft_x
        shifts_calc_y = rdisp[1]*0 + sft_y
        rot_angle = 0
        if sft_x and sft_y:
            vector_ratio = [np.median(shifts_destr_x/sft_x), np.median(shifts_destr_y/sft_y)]
        else:
            vector_ratio = [0.0, 0.0]
    
    disp_output = {'ref': dots_image_ref,
                   'sft': dots_image_sft,
                   'reg': data_sft_noisy_reg,
                   'rdisp': rdisp,
                   'disp': disp,
                   'rot_angle': rot_angle,
                   'input_shift': apply_shift,
                   'measured_shift': calc_shift,
                   'shifts_calc': [shifts_calc_x, shifts_calc_y],
                   'shifts_destr': [shifts_destr_x, shifts_destr_y]}
    return vector_ratio, disp_output
