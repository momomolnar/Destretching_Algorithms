
FUNCTION calc_synthetic_destretch, kernel_size, lowpass_factor=lowpass_factor, apod_percent=apod_percent, $
             input_image=input_image, output_image=output_image, regenerate_image=regenerate_image, dot_radius=dot_radius, dot_step=dot_step, $
             apply_rotation=apply_rotation, apply_shift=apply_shift, noise_level=noise_level, disp_output=disp_output, $
             debug_level=debug_level, plot_level=plot_level,flip_image=flip_image, gaussian_blur=gaussian_blur

IF N_ELEMENTS(dot_radius)       EQ 0 THEN dot_radius       = 5.1
IF N_ELEMENTS(dot_step)         EQ 0 THEN dot_step         = 12

IF N_ELEMENTS(regenerate_image) EQ 0 THEN regenerate_image = 0

IF N_ELEMENTS(apply_rotation)   EQ 0 THEN apply_rotation   = 0.0
IF N_ELEMENTS(apply_shift)      EQ 0 THEN apply_shift      = [0.0, 0.0]

IF N_ELEMENTS(noise_level)      EQ 0 THEN noise_level      = 0.0
IF N_ELEMENTS(debug_level)      EQ 0 THEN debug_level      = 0
IF N_ELEMENTS(plot_level)       EQ 0 THEN plot_level       = 0
IF N_ELEMENTS(flip_image)       EQ 0 THEN flip_image       = 0

IF N_ELEMENTS(gaussian_blur)    EQ 0 THEN gaussian_blur    = 1.0

;crosscor_peakfit_type = 'niblack_complex'
crosscor_peakfit_type = 'niblack_simple'

IF noise_level GT 0.0 THEN add_noise = 1 ELSE add_noise = 0

rot_angle    = apply_rotation
; this is the center of pixel of any rotation that is applied
; we are assuming we have a square image or at least the same value can be used for both axes
cent         = 500
; ratio of subfield spacing with respect to kernel size - passed through to the destretch procedure
spacrat      = 0.5

sft_x        = apply_shift[0]
sft_y        = apply_shift[1]

; this will be used in generating the synthetic image as a series of subfields
; the value should be bigger than the size of the radial distance subfields
; the radial distance subfields should be bigger than the maximum radius of the dot size 
buffer       = 25

; ----------------------------------------------------------------------
; if no input image is supplied, or new image generation is specifically requested, 
; then generate a new syntheic image
; ----------------------------------------------------------------------
IF N_ELEMENTS(input_image) LE 2 OR (regenerate_image EQ 1) THEN BEGIN
    numx      = 1000
    numy      = 1000
    nstep     = FIX(numx/dot_step)
    cent      = numx / 2.0
    ; magnitude of randomization - to scale the normally distributed random numbers
    ; the initial distance between dots is divided by this value, 
    ; so a bigger number means less randomness
    dot_loc_randomness = 2.0

    ; create an empty image array of desired size, plus buffer
    dots_image_ref = fltarr(numx + buffer*2, numy + buffer*2)                                                                                                                   
    dots_image_sft = fltarr(numx + buffer*2, numy + buffer*2)

    ; generate coordinates for dots, starting with regular grid
    ;xcoord = reform(rdisp[0,*,0])
    ;ycoord = reform(rdisp[1,0,*])
    xcoord = findgen(nstep) * numx/nstep
    ycoord = findgen(nstep) * numy/nstep
    xcoord += (numx - MAX(xcoord)) / 2.
    ycoord += (numy - MAX(ycoord)) / 2.
    ; now apply some randomization to coordinate locations
    seed0 = systime(/sec)                            
    seed0 = (seed0 * 1e5 - LONG64(seed0 * 1e5)) * 1e6
    coord_rand = randomn(seed0,[2, nstep, nstep]) * dot_step/dot_loc_randomness

    xcoord_num = N_ELEMENTS(xcoord)
    ycoord_num = N_ELEMENTS(ycoord)

    dot_coord_ref = FLTARR(2, xcoord_num, ycoord_num)

    ; sum regularly spaced dot postions with random jitter
    for xx=0, xcoord_num-1 do begin
        for yy=0, ycoord_num-1 do begin
            xcoord_jitter = (xcoord[xx] + coord_rand[0,xx,yy]) >0 <(numx-1)
            ycoord_jitter = (ycoord[yy] + coord_rand[1,xx,yy]) >0 <(numy-1)
            dot_coord_ref[*,xx,yy] = [xcoord_jitter, ycoord_jitter]
        endfor
    endfor

    ; if image rotation is being applied, calculate x-y shifts for each dot position
    IF apply_rotation THEN BEGIN
       beta               = atan(dot_coord_ref[0,*,*] - cent, dot_coord_ref[1,*,*] - cent)                                                                                  
        alpha              = rot_angle * !DTOR
        distance           = sqrt(TOTAL((dot_coord_ref - cent)^2,1))
        shifts_calc_x      = reform((sin(beta + alpha) * distance) - (sin(beta) * distance))
        shifts_calc_y      = reform((cos(beta + alpha) * distance) - (cos(beta) * distance))
    ENDIF

    ; compute dots for each position
    for xx=0, xcoord_num-1 do begin
        for yy=0, ycoord_num-1 do begin
             ; this was the original brute force approach, generating the radial distances over the full image for each dot
             ; that turned out to be slow... (and unnecessary) 
             ;dots_regular_1 += ABS(-(radial_distances([1,1000,1000], dot_coord_ref[*,xx,yy])<dot_radius - dot_radius))
             ;dots_regular_2 += ABS(-(radial_distances([1,1000,1000], dot_coord_ref[*,xx,yy] + [sft_x,sft_y])<dot_radius - dot_radius))

             ; select shifts to be applied to each dot position
             IF apply_rotation THEN BEGIN
                 offsets = [shifts_calc_x[xx,yy], shifts_calc_y[xx,yy]]
             ENDIF ELSE BEGIN
                 offsets = apply_shift
             ENDELSE

             ; now generate dots, initially as small sub-images, but then add those to the global image
             subfield_pix    = FIX(dot_coord_ref[*,xx,yy]) - buffer
             dot_subfield    = ABS(-(radial_distances([1, buffer*2, buffer*2], dot_coord_ref[*,xx,yy] - subfield_pix)<dot_radius - dot_radius))
             dots_image_ref[subfield_pix[0] + buffer:subfield_pix[0] + buffer*3-1, subfield_pix[1] + buffer:subfield_pix[1] + buffer*3-1] += dot_subfield

             dot_subfield    = ABS(-(radial_distances([1, buffer*2, buffer*2], dot_coord_ref[*,xx,yy] + offsets - subfield_pix)<dot_radius - dot_radius))
             dots_image_sft[subfield_pix[0] + buffer:subfield_pix[0] + buffer*3-1, subfield_pix[1] + buffer:subfield_pix[1] + buffer*3-1] += dot_subfield
        endfor
    endfor

    ; trim off buffer
    dots_image_ref = dots_image_ref[buffer:buffer+numx-1, buffer:buffer+numy-1]
    dots_image_sft = dots_image_sft[buffer:buffer+numx-1, buffer:buffer+numy-1]
    
ENDIF ELSE BEGIN
; ----------------------------------------------------------------------
; if an input image is supplied, then simply use cubic spline interpolation to introduce
; requested offsets
; ----------------------------------------------------------------------

    dots_image_ref = input_image
    scene_size     = size(dots_image_ref)
    numx           = scene_size[1]
    numy           = scene_size[2]

    IF apply_rotation THEN BEGIN
        dots_image_sft = ROT(dots_image_ref, rot_angle, 1.0, cent, cent, /pivot, cubic=-0.5)
    ENDIF ELSE BEGIN
        dots_image_sft = shift_im_cbc(dots_image_ref, apply_shift[0], apply_shift[1])
    ENDELSE
ENDELSE
; ----------------------------------------------------------------------

; the effective shift from the shift_im_cbc program might be different from the
; requested shift, so we will determine the actual shift using a cross-correlation
; over the whole image
calc_shift     = xyoff(dots_image_sft, dots_image_ref, 968, 968, /quiet)

; if requested, apply Gaussian smoothing
IF gaussian_blur GT 0 THEN BEGIN
    dots_image_ref = gauss_smooth(dots_image_ref, gaussian_blur, /edge_trunc)
    dots_image_sft = gauss_smooth(dots_image_sft, gaussian_blur, /edge_trunc)
ENDIF

if plot_level GE 1 then tvscl, dots_image_ref
;tvscl, data_ref - data_sft

; ----------------------------------------------------------------------
; add random (photon) noise if requested
IF add_noise EQ 1 THEN BEGIN
    seed1 = systime(/sec)                            
    seed2 = systime(/sec)                            
    seed1 = (seed1 * 1e5 - LONG64(seed1 * 1e5)) * 1e6
    seed2 = (seed2 * 1e5 - LONG64(seed2 * 1e5)) * 1e6

    data_ref_norm = dots_image_ref/MEAN(dots_image_ref)
    data_sft_norm = dots_image_sft/MEAN(dots_image_sft)

    dots_image_ref = (data_ref_norm) + (data_ref_norm * randomn(seed1,[numx, numy]) * noise_level)
    dots_image_sft = (data_sft_norm) + (data_sft_norm * randomn(seed2,[numx, numx]) * noise_level)
ENDIF
; ----------------------------------------------------------------------

; maybe the user want to flip the image in case to test whether there is some
; discrepancy in the amount of structure in the two array directions?
IF flip_image NE 0 THEN BEGIN
    dots_image_ref = ROTATE(dots_image_ref,4)
    dots_image_sft = ROTATE(dots_image_sft,4)
ENDIF

; if the user requested the output image, this will copy it into their supplied array
output_image = dots_image_ref

; ----------------------------------------------------------------------
; now run destretching on reference and shifted image - finally!
; ----------------------------------------------------------------------

data_sft_noisy_reg = destretch(dots_image_sft, dots_image_ref, bytarr(kernel_size, kernel_size),rdisp=rdisp, disp=disp, $
                spacing_ratio=spacrat, crosscor_fit= crosscor_peakfit_type, debug_level=debug_level, lowpass_factor=lowpass_factor, $
                apod_percent=apod_percent, plot_level=plot_level)

shifts_destr_x = (disp[0,*,*] - rdisp[0,*,*])
shifts_destr_y = (disp[1,*,*] - rdisp[1,*,*])

; ----------------------------------------------------------------------

; a final summary calculation comparing the ratio of expected (input) displacements to those actually measured
IF apply_rotation THEN BEGIN
    ; recalculate rotation-based shifts
    beta               = atan(rdisp[0,*,*] - cent, rdisp[1,*,*] - cent)                                                                                  
    alpha              = rot_angle * !DTOR
    distance           = sqrt(TOTAL((rdisp - cent)^2,1))
    shifts_calc_x      = (sin(beta + alpha) * distance) - (sin(beta) * distance)
    shifts_calc_y      = (cos(beta + alpha) * distance) - (cos(beta) * distance)

    vector_ratio = [ median(shifts_destr_x / shifts_calc_x) , $
      median(shifts_destr_y / shifts_calc_y) ]
ENDIF ELSE BEGIN
    ; use the rdisp array to get an array of the proper size, but then zero it out and add back in the input shifts
    shifts_calc_x      = rdisp[0,*,*] * 0.0 + sft_x
    shifts_calc_y      = rdisp[1,*,*] * 0.0 + sft_y
    rot_angle          = 0

    IF (sft_x NE 0) AND (sft_y NE 0) THEN $ 
        vector_ratio = [ median(shifts_destr_x / sft_x) , median(shifts_destr_y / sft_y) ] $
        ELSE vector_ratio = [0.0,0.0]
ENDELSE

; ----------------------------------------------------------------------

; generate structure with destretch parameters
disp_output        = create_struct('rdisp', rdisp, $
                                   'disp',  disp, $
                                   'rot_angle',  rot_angle, $
                                   'input_shift', apply_shift, $
                                   'measured_shift', calc_shift, $
                                   'shifts_calc',  [shifts_calc_x, shifts_calc_y], $
                                   'shifts_destr', [shifts_destr_x, shifts_destr_y])


RETURN,vector_ratio 

end