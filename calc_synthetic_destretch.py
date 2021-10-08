#!/usr/bin/env python3

from destretch import *
import time
import tqdm
import matplotlib.pyplot as plt
from astropy.io import fits
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
    
    if (not np.any(input_image)) or regenerate_image:
        if debug_level:
            print('Generating random images')
        # Generate image of random dots
        numx = 1000
        numy = 1000
        nstep = int(numx / dot_step)
        dot_loc_randomness = 2.0
        dots_image_ref = np.zeros((numy + buffer*2, numx + buffer*2))
        dots_image_sft = np.zeros((numy + buffer*2, numx + buffer*2))
        xcoord = np.arange(nstep) * numx/nstep
        ycoord = np.arange(nstep) * numy/nstep
        xcoord += (numx - xcoord.max()) / 2.
        ycoord += (numy - ycoord.max()) / 2.
        rng = np.random.default_rng()
        coord_rand = rng.normal(size=(2, nstep, nstep)) * dot_step/dot_loc_randomness
        coord_rand = np.array(coord_rand)
        xcoord_num = xcoord.size
        ycoord_num = ycoord.size
        dot_coord_ref = np.zeros((2, ycoord_num, xcoord_num))
        
        for y in range(ycoord_num):
            for x in range(xcoord_num):
                xcoord_jitter = np.clip(xcoord[x] + coord_rand[0, y, x], 0, numx-1)
                ycoord_jitter = np.clip(ycoord[y] + coord_rand[1, y, x], 0, numy-1)
                dot_coord_ref[:, y, x] = [ycoord_jitter, xcoord_jitter]
        
        if apply_rotation:
            if debug_level:
                print('Rotating random images')
            beta = np.arctan2(dot_coord_ref[0] - cent, dot_coord_ref[1] - cent)
            alpha = rot_angle * np.pi/180.
            distance = np.sqrt(np.sum((dot_coord_ref - cent)**2, 0))
            shifts_calc_x = np.sin(beta + alpha) * distance - np.sin(beta) * distance
            shifts_calc_y = np.cos(beta + alpha) * distance - np.cos(beta) * distance
        
        for y in range(ycoord_num):
            for x in range(xcoord_num):
                if apply_rotation:
                    offsets = [shifts_calc_x[y, x], shifts_calc_y[y, x]]
                else:
                    offsets = apply_shift
                subfield_pix = np.around(dot_coord_ref[:, y, x]).astype(int) - buffer
                dot_subfield = np.clip(dot_radius - radial_distances([buffer*2], dot_coord_ref[:, y, x] - subfield_pix), 0, None)
                dots_image_ref[subfield_pix[1] + buffer:subfield_pix[1] + buffer*3, subfield_pix[0] + buffer:subfield_pix[0] + buffer*3] += dot_subfield
                dot_subfield = np.clip(dot_radius - radial_distances([buffer*2], dot_coord_ref[:, y, x] + offsets - subfield_pix), 0, None)
                dots_image_sft[subfield_pix[1] + buffer:subfield_pix[1] + buffer*3, subfield_pix[0] + buffer:subfield_pix[0] + buffer*3] += dot_subfield
        dots_image_ref = dots_image_ref[buffer:buffer+numy, buffer:buffer+numx]
        dots_image_sft = dots_image_sft[buffer:buffer+numy, buffer:buffer+numx]
        if debug_level:
            print('Images successfully created')
    else:
        if debug_level:
            print('Using provided images')
        dots_image_ref = input_image
        scene_size = dots_image_ref.shape
        numy = scene_size[0]
        numx = scene_size[1]
        
        if apply_rotation:
            dots_image_sft = rotate_around_point(dots_image_ref, rot_angle, [cent, cent], 3)
        else:
            dots_image_sft = shift(dots_image_ref, apply_shift)
    
#    corr = signal.correlate2d(dots_image_sft, dots_image_ref, boundary='symm')
    fs = np.fft.fft2(dots_image_sft)
    fr = np.fft.fft2(dots_image_ref)
    f = np.abs(np.fft.fftshift(np.fft.ifft2(fr.conj()*fs)))
    fx = np.max(f)
    ip = np.where(f==fx)
    iy = ip[0]
    ix = ip[1]
    pp = np.squeeze(np.array(ip)[::-1].astype(float))
    if ix and (numx - ix - 1):
        pp[0] = pp[0] + 0.5*(f[iy,ix+1] - f[iy,ix-1])/(2*fx - f[iy,ix-1] - f[iy,ix+1])
    if iy and (numy - iy - 1):
        pp[1] = pp[1] + 0.5*(f[iy+1,ix] - f[iy-1,ix])/(2*fx - f[iy-1,ix] - f[iy+1,ix])
    calc_shift = pp - 500
    if debug_level:
        print('input_shift = ' + repr(apply_shift))
        print('calc_shift = ' + repr(calc_shift))
    if gaussian_blur:
        if debug_level:
            print('Degrading with gaussian blur')
        dots_image_ref = gaussian_filter(dots_image_ref, gaussian_blur)
        dots_image_sft = gaussian_filter(dots_image_sft, gaussian_blur)
    else:
        if debug_level:
            print('No blurring')
    if plot_level:
        if debug_level:
            print('Plotting reference image')
        plt.ion()
        plt.imshow(dots_image_ref)
    else:
        if debug_level:
            print('No plotting')
    
    if add_noise:
        if debug_level:
            print('Degrading image with shot noise')
        rng2 = np.random.default_rng()
        rng3 = np.random.default_rng()
        data_ref_norm = dots_image_ref / np.mean(dots_image_ref)
        data_sft_norm = dots_image_sft / np.mean(dots_image_sft)
        dots_image_ref = data_ref_norm * (1 + rng2.normal(size=(numy, numx))*noise_level)
        dots_image_sft = data_sft_norm * (1 + rng3.normal(size=(numy, numx))*noise_level)
    else:
        if debug_level:
            print('No shot noise')
    
    if flip_image:
        if debug_level:
            print('Transposing images')
        dots_image_ref = dots_image_ref.T
        dots_image_sft = dots_image_sft.T
    output_image = dots_image_ref
    
    if debug_level:
        print('Now destretching shifted image, please wait...')
    data_sft_noisy_reg, disp, rdisp, d_info = reg(dots_image_sft,
                                                  dots_image_ref,
                                                  kernel_size,
                                                  mf=apod_percent)
    if debug_level:
        print('Destretching complete!')
    shifts_destr_x = (disp - rdisp)[0]
    shifts_destr_y = (disp - rdisp)[1]
    
    if apply_rotation:
        beta = np.arctan2(rdisp[1] - cent, rdisp[0] - cent)
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
    
    
if __name__ == "__main__":

    sft_range = 4

    wl_mfbd_in = fits.open('/home/ryan/Destretching_Algorithms/test/whitelight.burst.20170423_172226.seq00.im000.rot.fits')[0].data.astype(float)
    wl_mfbd = np.zeros((1000, 1000)) + wl_mfbd_in.mean()
    wl_mfbd[8:-8, 8:-8] = wl_mfbd_in[8:-8, 8:-8]
    wl_mfbd /= wl_mfbd.mean()
    
    num_repeats = 100
    
    rng = np.random.default_rng()
    shifts = rng.random(size=(2, num_repeats)) * sft_range * 2 - sft_range
    
    kerns = [16,24,32,48,64,96,128,192]
    apods = (np.arange(25) + 1)*0.01
    dotradii = np.arange(15) + 1
    smoothing = [0.1,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0]
    
    vary_kernels = 1
    vary_apod    = 0
    vary_dotrad  = 0
    vary_smooth  = 0

    use_real_im = 0
    
    if vary_kernels:
        params_num = len(kerns)
        indep_var = kerns
    elif vary_apod:
        params_num = len(apods)
        indep_var = apods
    elif vary_dotrad:
        params_num = len(dotradii)
        indep_var = dotradii
    elif vary_smooth:
        params_num = len(smoothing)
        indep_var = smoothing
    
    destr_stats_byparam = np.zeros((19, params_num))
    disp_stats_sfts_all = np.zeros((10, num_repeats, params_num))
    
    for param_idx in range(params_num):
        kern_size = 32
        dotrad = 6
        apod_val = 0.1
        smooth_val = 1.0
        
        if vary_kernels:
            kern_size = kerns[param_idx]
        elif vary_apod:
            apod_val = apods[param_idx]
        elif vary_dotrad:
            dotrad = dotradii[param_idx]
        elif vary_smooth:
            smooth_val = smoothing[param_idx]
        
        if use_real_im:
            input_image = wl_mfbd
        else:
            input_image = 0
        
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
              "  ; kernel = ", kern_size,
              " ; apodization = ", apod_val,
              " ; dot radius = ", dotrad,
              " ; gaussian blurring = ", smooth_val)
        start_time = time.time()
        
        disp_stats_shifts = np.zeros((10, num_repeats))
        
        for nn in tqdm.tqdm(range(num_repeats)):
            sfts = shifts[:, nn]
            
            debug_level = 1 if nn == 0 else 0
            plot_level = 2 if nn == 0 else 0
            
            gg2, disp_out_temp = calc_synthetic_destretch(kern_size,
                                                          apod_percent=apod_val,
                                                          apply_shift=sfts,
                                                          dot_radius=dotrad,
                                                          dot_step=7,
                                                          input_image=input_image,
                                                          debug_level=debug_level,
                                                          plot_level=plot_level,
                                                          flip_image=0,
                                                          gaussian_blur=smooth_val)
            
            stats = [disp_out_temp['measured_shift'][0],
                     disp_out_temp['measured_shift'][1],
                     np.mean(disp_out_temp['shifts_destr'][0]),
                     np.mean(disp_out_temp['shifts_destr'][1]),
                     np.median(disp_out_temp['shifts_destr'][0]),
                     np.median(disp_out_temp['shifts_destr'][1]),
                     np.std(disp_out_temp['shifts_destr'][0]),
                     np.std(disp_out_temp['shifts_destr'][1]),
                     np.mean(disp_out_temp['ref']),
                     np.std(disp_out_temp['ref'])]
            disp_stats_shifts[:, nn] = stats
        disp_stats_sfts_all[:, :, param_idx] = disp_stats_shifts
#        print(repr(disp_stats_shifts))
        
        linfit_sfts_mean_x = np.polyfit(disp_stats_shifts[0, :], disp_stats_shifts[2, :], 1)
        linfit_sfts_mean_y = np.polyfit(disp_stats_shifts[1, :], disp_stats_shifts[3, :], 1)
        linfit_sfts_med_x = np.polyfit(disp_stats_shifts[0, :], disp_stats_shifts[4, :], 1)
        linfit_sfts_med_y = np.polyfit(disp_stats_shifts[1, :], disp_stats_shifts[5, :], 1)

        sfts_rms_med_x = np.median(disp_stats_shifts[6, :])
        sfts_rms_med_y = np.median(disp_stats_shifts[7, :])
        
        linfit_rms_x = np.polyfit(np.abs(disp_stats_shifts[0, :]), disp_stats_shifts[6, :], 1)
        linfit_rms_y = np.polyfit(np.abs(disp_stats_shifts[1, :]), disp_stats_shifts[7, :], 1)

        image_int_ave      = np.median(disp_stats_shifts[8, :])
        image_rms_ave      = np.median(disp_stats_shifts[9, :])

        image_sft_change_x = np.median(disp_stats_shifts[0, :] / shifts[0, :])
        image_sft_change_y = np.median(disp_stats_shifts[1, :] / shifts[1, :])

        end_time = time.time()
        time_per_cycle = (end_time - start_time) / num_repeats
        
        destr_stats_byparam[:, param_idx] = [linfit_sfts_mean_x[1],
                                             linfit_sfts_mean_x[0],
                                             linfit_sfts_mean_y[1],
                                             linfit_sfts_mean_y[0],
                                             linfit_sfts_med_x[1],
                                             linfit_sfts_med_x[0],
                                             linfit_sfts_med_y[1],
                                             linfit_sfts_med_y[0],
                                             sfts_rms_med_x,
                                             sfts_rms_med_y,
                                             linfit_rms_x[1],
                                             linfit_rms_x[0],
                                             linfit_rms_y[1],
                                             linfit_rms_y[0],
                                             image_int_ave,
                                             image_rms_ave,
                                             time_per_cycle,
                                             image_sft_change_x,
                                             image_sft_change_y]
        print(repr(destr_stats_byparam[:, param_idx]))
    
    for param in range(19):
        print(repr(destr_stats_byparam[param]))
    
    fig, ax = plt.subplots(5, 4, figsize=(16, 10))
    for i in range(5):
        for j in range(4):
            try:
                dep_var = np.abs(destr_stats_byparam[4*i+j])
                if np.median(dep_var) > 0.8 and np.max(dep_var) < 1.2:
                    dep_var = max(1, np.max(dep_var)*1.0001) - dep_var
                ax[i, j].scatter(indep_var, dep_var)
                ax[i, j].set_xscale('log')
                ax[i, j].set_yscale('log')
            except:
                break
    plt.tight_layout()
    plt.savefig('synthetic_destretch_results_{}{}{}{}.pdf'.format(vary_kernels,
                                                                  vary_apod,
                                                                  vary_dotrad,
                                                                  vary_smooth))
