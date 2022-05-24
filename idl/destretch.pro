;	reg.cpp		register a sequence of 1 or more images
;			(uses fft's to compute 2d cross corrrelations)
;	22 Jun 92 phw&tr wx,wy always even
;
;	22 Jun 92 phw&tr Include X as a workstation device;
;			save !x,!y sys variables in common, not scratch file
;	04 Feb 91 phw	correct a bug in repair (kkx, instead of kx, etc);
;			correct a buglet in mkcps (debug might be undef'd)
;	18 Jan 91 phw	add repair of control point displacements;
;			add debug to reg_com
;	28 Nov 90 phw	(bilin) speedup;
;			(bspline) better limit checks on t0,tn,s0,sn.
;	20 Nov 90 phw	(bspline) messing with boundry, so it shifts;
;			(cps) correct undef'd rcps.
;	19 Nov 90 phw	(bspline) fix boundry error introduced 25 Sep
;	26 Oct 90 phw	add setup & undo;
;			save & restore system variables
;	27 Sep 90 phw	change cpstruct to common reg_comm;
;			use rotate (x,2) instead of reverse (x);
;			reorder so only one loading required.
;	26 Sep 90 phw	fine tune mkcps & bspline
;	25 Sep 90 phw	kernal wanders in larger area; wx, wy
;	06 Apr 90 phw	usage comments added;
;			irrelevant interlace stuff deleted.
;	19 sep 89 phw	rename 'decide' to 'mkcps'
;			make callable without cpstruct
;	18 sep 89 phw	generalize mask for unequal kx,ky
;	06 sep 89 phw	include B-spline, bilin
;	29 aug 89 phw	fix center points; trying B-splines
;	28 aug 89 phw	masking of reference
;	15 aug 89 phw	remove qwx stuff; add rms of fit
;	10 aug 89 phw	modified for fft method
;	early 89 phw	written

;	This package of IDL procedures and functions demonstrates
; one way that a sequence of images can be registered and destretched.
; The method used is to tile the image(s) with (usually square) kernels.
; For each kernel:
;	a 2D cross correlation is computed with respect to a reference image;
;	the position of the maximum is interpolated for sub-pixel position;
;	the positions are collected into a table of control (tile) points;
;	a B-spline surface is fitted to these tile points; and finally
;	the scene is resampled on the non-uniform grid.
;
;	The above may be accomplished with the function 'reg' (below).
; For example, the following IDL code should do the job:
;
;	...
;	sequence = bytarr (256,256,200)	; series of 2D arrays to destretch
;	readu, unit, sequence		; get'um into memory
;	ref = sequence (*,*,0)		; select reference frame
;	kernel = bytarr (15,15)		; select kernel (tile) size
;	destretched = reg (sequence, ref, kernel) ; doit - result will be
;	...				; same size fltarr as input sequence
;
;	We have used the method iteratively, with successively
; smaller kernel sizes (suggest 64^2, 40^2, 25^2 for 256^2 images).  We
; have also had better luck if the entire image is first correlation
; tracked.  More recently (26 Sep 90) refinements have been added
; that permit successful use of a small kernal, one time (i.e.,
; iteration should be unnecessary).
;
;	The implementation consists of several routines that are intended
; to be user callable, to give the user control over such things as
; updateing the reference scene.  For this, the reader's attention is
; called to the routines:
;
;	mkcps	to compute reference control point coordinates (one time)
; 	cps	to compute the offsets for current scene from reference
;	repair	to repair control point displacements
;	doreg	to apply destretch to scene, which may become the
;		new reference.
;
; 	Debug plots are generated and various debug print statements
; to trace the progress through the package.  The package works, but
; somewhat ungracefully on a tektronix compatible terminal, as well
; as on a sun workstation. 
; The plots are not crucial, though they give some indication of how
; well the tracking is doing.
;
; For more information contact:
;
;	Phil Wiborg
;	National Solar Observatory
;	Sun Spot, NM 88349
;	505-434-7000
;	pwiborg@sunspot.noao.edu
; or
;	Thomas Rimmele
;	Kiepenheuer Institut f"ur 
;	Sonnenphysik
;	7800 Freiburg
;	0761-3198-0
;	tr@kis.uni-freiburg.de
;
; **********************************************************
; ********************  PRO: define_destr_info  *******************
; **********************************************************;+
; PURPOSE:
;  This function generates a blank destr_info array
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;  rcps_size = size of control_point array to include in structure
; OUTPUTS:
;  destr_info: shared variable containing destretching parameters and control points
;
; MODIFICATION HISTORY:
;  August 2021: K. Reardon - cleaned up and modernized
;-
FUNCTION define_destr_info, rcps_size=rcps_size

IF N_ELEMENTS(rcps_size) GE 1 THEN $
     rcps_out = FLTARR(rcps_size[0], rcps_size[1], rcps_size[2]) $
ELSE rcps_out = 0

destr_info = Create_Struct(    'kx',          0, $
                               'ky',          0, $
                               'spacing_x',   0, $
                               'spacing_y',   0, $
                               'border_x',    0, $
                               'border_y',    0, $
                               'wx',          0, $
                               'wy',          0, $
                               'bx',          0, $
                               'by',          0, $
                               'cpx',         0, $
                               'cpy',         0, $
                               'ref_sz_x',    0, $
                               'ref_sz_y',    0, $
                               'scene_sz_x',  0, $
                               'scene_sz_y',  0, $
                               'scene_sz_z',  0, $
                               'rcps',        rcps_out, $
                               'subfield_correction', '', $
                               'max_fit_method',      '', $
                               'do_plots',    0, $
                               'debug',       0)

return, destr_info

end

; **********************************************************
; ******************** FUNCTION: bilin  *******************
; **********************************************************
;+
; PURPOSE:
;  generate a bilinear interpolation of a scene
;
; INPUTS:
;  scene: input image
;
; KEYWORD PARAMETERS:
; 
; OUTPUTS:
;
; MODIFICATION HISTORY:
;  26 October 1990: Originally written by Phil Wiborg
;  August 2021: K. Reardon - cleaned up and modernized
;-
FUNCTION	bilin, scene, coords_new, bilinear=bilinear
;
;	bilinear interpolation of scene, s
;

IF N_ELEMENTS(bilinear) EQ 0 THEN bilinear=1

if NOT bilinear then begin	
   ; this is a simple nearest neighbor

    x = fix(coords_new(*,*,0) + 0.5)
    y = fix(coords_new(*,*,1) + 0.5)
    scene_destr = scene(x,y)

endif else begin		
    ; bilinear

    x = coords_new(*,*,0)
    y = coords_new(*,*,1)

    x0 = fix(x)
    x1 = fix(x+1)
    y0 = fix(y)
    y1 = fix(y+1)

    ; find the fractional part of the modified coordinates
    fx = x mod 1.
    fy = y mod 1.

    scene_flt = float(scene)

; original (slow) version
;    ssfx = (scene_flt(x1,y0)-scene_flt(x0,y0))*fx
;    scene_destr = ssfx + $
;	(scene_flt(x0,y1) - scene_flt(x0,y0))*fy + $
;	((scene_flt(x1,y1) - scene_flt(x0,y1))*fx - ssfx)*fy + $
;	scene_flt(x0,y0)

; optimized version
    ss00 = scene_flt(x0,y0)
    ss01 = scene_flt(x0,y1)
    ssfx00 = (scene_flt(x1,y0) - ss00)*fx
    ssfx01 = (scene_flt(x1,y1) - ss01)*fx
    scene_destr  = ss00 + ssfx00 + (ss01 - ss00 + ssfx01*fx - ssfx00)*fy
    endelse

return, scene_destr
end

;
; **********************************************************
; ******************** FUNCTION: extend  *******************
; **********************************************************
; PURPOSE:
;     extend reference and actual displacements to cover area beyond scene;
; INPUTS:
;   rdisp: fltarr(cpx,cpy) -  reference control points
;   disp:  fltarr(cpx,cpy) -  actual displaced position of control points
;
; KEYWORD PARAMETERS:
;  extend_with_offsets: extend displacement array by replicating measured displacements
;                       at borders of measured array. Set to 0 to set all extended
;                       displacements to 0
; OUTPUTS:
;   rdisp_ext   = array of reference control points, extended beyond measured range
;   offsets_ext = array of displacements, extended beyond measured range
; NOTES:
;
; MODIFICATION HISTORY:
;  26 October 1990: Originally written by Phil Wiborg
;  August 2021: K. Reardon - cleaned up and modernized
;-
pro	extend, rdisp_1d, disp_1d, rdisp_ext, offsets_ext, extend_with_offsets= extend_with_offsets	

IF N_ELEMENTS(extend_with_offsets) EQ 0 THEN extend_with_offsets=1

num_extpts = 3
; get size of displacement array
;   add num_extpts control points on each side
disp_sz = size(disp_1d)
disp_nx = disp_sz(1) + num_extpts*2
disp_ny = disp_sz(2) + num_extpts*2

rdisp_ext = fltarr(disp_nx, disp_ny)
step_x    = rdisp_1d(1,0) - rdisp_1d(0,0)
; define the starting point, which is num_extpts step sizes before the input position
start_pt  = rdisp_1d(0,0) - num_extpts * step_x
; generate a new array of evenly spaced control points
new_steps = findgen(disp_nx) * step_x + start_pt
; and fill the new rdisp array with those control points
for j=0, disp_ny-1 do rdisp_ext(*,j) = new_steps

IF extend_with_offsets THEN BEGIN
    ; extend displacement array by replicating the values of the displacements along the 
    ;   borders of measured array. This ensures some consistencies in the values

    ; create an extended array and populate the center portion
    ; with the measured offsets
    offsets_ext    = fltarr (disp_nx, disp_ny)
    offsets_ext[num_extpts:disp_nx-num_extpts-1, num_extpts:disp_ny-num_extpts-1] = disp_1d - rdisp_1d

    ; take the bottom row of measured offset and replicate it into the extended array of offsets below
    ;   note: this could all be done with rebin instead of for loops...
    x = offsets_ext(*, num_extpts)
    for extra_row=0, num_extpts-1 do offsets_ext(*, extra_row) = x

    ; take the top row of measured offset and replicate it into the extended array of offsets above
    x = offsets_ext(*, disp_ny-num_extpts-1)
    for extra_row=disp_ny-num_extpts,disp_ny-1 do offsets_ext(*, extra_row) = x

    ;rotate array and repeat filling process
    offsets_ext = transpose (offsets_ext)

    ; take the bottom row of measured offset and replicate it into the extended array of offsets below
    x = offsets_ext(*, num_extpts)
    for extra_row=0, num_extpts-1 do offsets_ext(*, extra_row) = x

    ; take the top row of measured offset and replicate it into the extended array of offsets above
    x = offsets_ext(*, disp_ny-num_extpts-1)
    for extra_row=disp_ny-num_extpts,disp_ny-1 do offsets_ext(*, extra_row) = x

    ;rotate array back to original orientation
    offsets_ext = transpose (offsets_ext)

    ; add the extended array of offsets back to the extened array of reference control point positions
    offsets_ext = offsets_ext + rdisp_ext
ENDIF ELSE BEGIN
    ; alternatively, set all extended displacements to 0, by simply setting all the displacement values
    ;   in the extended portion of the array to the same values as reference positions
    offsets_ext = rdisp_ext
    offsets_ext[num_extpts:disp_nx-num_extpts-1, num_extpts:disp_ny-num_extpts-1] = disp_1d 
ENDELSE

return

end
; **********************************************************
; ******************** FUNCTION: bspline  ******************
; **********************************************************
; PURPOSE:
;  destretch scene using B-splines
;
; INPUTS:
;   scene: fltarr(nx,ny)     -  image to be destretched
;   rdisp: fltarr(2,cpx,cpy) -  reference control points
;   disp:  fltarr(2,cpx,cpy) -  actual displaced position of control points
;
; KEYWORD PARAMETERS:
; OUTPUTS:
;   coords_destr =  coordinates for destretched scene
; NOTES:
;     based on Foley & Van Dam: pp 521-536.
;
; MODIFICATION HISTORY:
;  26 October 1990: Originally written by Phil Wiborg
;  August 2021: K. Reardon - cleaned up and modernized
;  May 2022: K. Reardon - changed to a standalone function, as one interpolation
;                         option to be called by rescale_distorion_map
;-
FUNCTION bspline, rdisp, disp, destr_info=destr_info, do_plots=do_plots, $
                  extend_with_offsets=extend_with_offsets

IF N_ELEMENTS(do_plots) EQ 0 THEN do_plots=0   
IF N_ELEMENTS(extend_with_offsets) EQ 0 THEN extend_with_offsets=1

; this tells us the target size of the output array
scene_x = destr_info.scene_sz_x
scene_y = destr_info.scene_sz_y
ans     = fltarr (scene_x, scene_y, 2)

grid_spacing_x   = rdisp(0,1,0)-rdisp(0,0,0)
grid_spacing_y   = rdisp(1,0,1)-rdisp(1,0,0)

disp_sz = size(disp)
disp_nx = disp_sz(2)
disp_ny = disp_sz(3)

grid_extend_factor = 3
; find the biggest distance from the outside control point to the edge of the full scene
border_distance_x = max([min(rdisp[0,*,*],max=maxval_x),scene_x - maxval_x])
border_distance_y = max([min(rdisp[1,*,*],max=maxval_y),scene_y - maxval_y])
; define the number of extra "control points" that are needed on each edge
; - the number of extra pixels needed, divided by the size of the grid spacing
extend_range_x    = FIX(border_distance_x * grid_extend_factor / grid_spacing_x) + 1
extend_range_y    = FIX(border_distance_y * grid_extend_factor / grid_spacing_y) + 1

extend, reform(rdisp[0,*,*], disp_nx, disp_ny), $
            reform(disp[0,*,*], disp_nx, disp_ny), Rx, Px, $
            extend_with_offsets=extend_with_offsets, num_extpts=extend_range_x
extend, transpose(reform(rdisp[1,*,*], disp_nx, disp_ny)), $
            transpose(reform(disp[1,*,*], disp_nx, disp_ny)), Ry, Py, $
            extend_with_offsets=extend_with_offsets, num_extpts=extend_range_y
Ry = transpose (Ry)
Py = transpose (Py)

rdisp_extend_nx = n_elements(Rx[*,0])
rdisp_extend_ny = n_elements(Ry[0,*])

Ms  = [-1,3,-3,1, 3,-6,0,4, -3,3,3,1, 1,0,0,0]/6.
Ms  = reform (Ms, 4,4)
MsT = transpose(Ms)

; step through each control point in the extended array of reference points,
;     first in the y-axis
for v=0, rdisp_extend_ny-3 do begin
    t0 = Ry(1,v+1)
    tn = Ry(1,v+2)
    ; if the value of the pixel position at the current control point is outside
    ;     the bounds of the full array, we don't need to interpolate a value 
    ;     at those points (since they are outside the bounds of the final image)
    if (tn le 0) or (t0 ge scene_y-1) then begin
        ; nothing to do, outside of the final array, skip this iteration
    endif else begin
        t0 = max ([t0, 0])
        tn = min ([tn, scene_y])
        t = findgen(FIX(tn)-FIX(t0))/grid_spacing_y + (t0-Ry(1,v+1))/grid_spacing_y

        for u=0, rdisp_extend_nx-3 do begin
            s0 = Rx(u+1,v+1)
            sn = Rx(u+2,v+1)
            if (sn le 0) or (s0 ge scene_x-1) then begin
                ; nothing to do, skip this iteration
            endif else begin
                s0 = max ([s0,0])
                sn = min ([sn, scene_x])
                s = findgen(FIX(sn)-FIX(s0))/grid_spacing_x + (s0-Rx(u+1,v+1))/grid_spacing_x
                ;help,Ms,Px(u:u+3,v:v+3),MsT
                compx = reform (Ms # Px(u:u+3,v:v+3) # MsT, 4, 4)
                compy = reform (Ms # Py(u:u+3,v:v+3) # MsT, 4, 4)
                ; previously the secondary "patch" procedure was called to calculate
                ; the interpolated shifts
                ;ans(s0:sn-1,t0:tn-1,*) = patch (compx, compy, s, t)
                
                ; now we just move that algorithm here
                ss = reform([s^3,s^2,s,replicate(1.,n_elements(s))],n_elements(s),4)
                tt = transpose(reform([t^3,t^2,t,replicate(1.,n_elements(t))],n_elements(t),4))
                a1x = ss # compx
                a1y = ss # compy
                a2x = a1x # tt
                a2y = a1y # tt

                ans(s0:sn-1,t0:tn-1,0) = ss # compx # tt
                ans(s0:sn-1,t0:tn-1,1) = ss # compy # tt
       
                if (do_plots GE 2) and (u gt 0) and (v gt 0) then begin
                    plots, reform(ans(s0,t0,*),2,1), /dev, psym=7
                    wait,0.
                endif

        endelse    ; nextu:
        endfor

    endelse  ; nextv:
endfor

return, ans
end

; ************************************************************************
; ******************** FUNCTION: rescale_distorion_map  ******************
; ************************************************************************
; PURPOSE:
;  interpolate shifts measured at each control point to a grid of 
;  distortions at each pixel in the full scene
;
; INPUTS:
;    rdisp: fltarr(2,cpx,cpy) -  reference control points
;    disp:  fltarr(2,cpx,cpy) -  actual displaced position of control points
;
; KEYWORD PARAMETERS:
;    destr_info: structure containing information on destretch parameters
;    do_plots: plot some outputs of the interpolation
;    original_version: use the original B-spline algorithm; default is to use
;                      IDL's built in INTERPOLATE command
;    scale_fact: scale factor to change the amplitude of the measured distortions. 
;                probably needed previously due to deficiencies in interpolation
;                or cross-correlation. Default of 1.0 is probably fine now.
; OUTPUTS:
;    coords_destr =  coordinates for destretched scene
; NOTES:
;
; MODIFICATION HISTORY:
;  26 October 1990: Originally written by Phil Wiborg
;  August 2021: K. Reardon - cleaned up and modernized
;  May 2022: K. Reardon - added INTERPOLATE option since it better maintains
;                         the small scale features of the measured distortions
;-
FUNCTION rescale_distorion_map, rdisp, disp, destr_info=destr_info, do_plots=do_plots, $
                                original_version=original_version, scale_fact=scale_fact

IF N_ELEMENTS(do_plots) EQ 0         THEN do_plots=0   
IF N_ELEMENTS(original_version) EQ 0 THEN original_version=0   
IF N_ELEMENTS(scale_fact) EQ 0       THEN scale_fact = 1.0

; amplify (or reduce) the size of the distortions to be remapped
disp_scaled = (disp - rdisp)*scale_fact + rdisp

grid_spacing_x = rdisp(0,1,0)-rdisp(0,0,0)
grid_spacing_y = rdisp(1,0,1)-rdisp(1,0,0)

disp_sz = size (disp)
disp_nx = disp_sz(2)
disp_ny = disp_sz(3)
    
scene_sz = size(scene)
scene_x = destr_info.scene_sz_x
scene_y = destr_info.scene_sz_y

; distortion map needs to be extended at borders so interpolation scheme can
; extrapolate distortions past the boundaries of the interior control points
; where the distortions were measured. This is done calling the "extend" function below.
;
; there are two ways to define the distortions extended at the edges:
; 1) one is to extend the reference positions out to new points. This essentially 
;       sets any subfield shifts around the boundaries to zero
;extend_with_offsets = 0       ; exterior control points set to zero displacements
; 2) make the distortions at the boundaries equal to those at the nearby 
;       control points, avoiding any discontinuities in the offsets
extend_with_offsets = 1        ; exterior control points drift with interior (best)

ans = fltarr (scene_x, scene_y, 2)

if original_version EQ 1 then begin
    if (destr_info.debug GE 1) then print,"Scaling shifts using orginal B-spline"

    ; splitting all the old-school interpolation out into a separate routine
    ans = bspline(rdisp, disp, destr_info=destr_info, do_plots=do_plots, $
                      extend_with_offsets=extend_with_offsets)
endif else begin
    if (destr_info.debug GE 1) then print,"Scaling shifts using IDL INTERPOLATE"

    grid_extend_factor = 3
    ; find the biggest distance from the outside control point to the edge of the full scene
    border_distance_x = max([min(rdisp[0,*,*],max=maxval_x),scene_x - maxval_x])
    border_distance_y = max([min(rdisp[1,*,*],max=maxval_y),scene_y - maxval_y])
    ; define the number of extra "control points" that are needed on each edge
    ; - the number of extra pixels needed, divided by the size of the grid spacing
    extend_range_x    = FIX(border_distance_x * grid_extend_factor / grid_spacing_x) + 1
    extend_range_y    = FIX(border_distance_y * grid_extend_factor / grid_spacing_y) + 1
    
    extend, reform(rdisp[0,*,*], disp_nx, disp_ny), $
            reform(disp_scaled[0,*,*], disp_nx, disp_ny), Rx, Px, $
            extend_with_offsets=extend_with_offsets, num_extpts=extend_range_x
    ; use transpose to make y-direction arrays comparable to x-direction arrays
    extend, transpose(reform(rdisp[1,*,*], disp_nx, disp_ny)), $
            transpose(reform(disp_scaled[1,*,*], disp_nx, disp_ny)), Ry, Py, $
            extend_with_offsets=extend_with_offsets, num_extpts=extend_range_y
    ; then transpose them back
    Ry = transpose (Ry)
    Py = transpose (Py)

    grid_start_x              = rdisp[0,0,0]
    grid_start_y              = rdisp[1,0,0]
    ;grid_spacing_x            = rdisp_x[1,0] - rdisp_x[0,0]
    ;grid_spacing_y            = rdisp_y[0,1] - rdisp_y[0,0]
    
    ; this defines a grid of points for which we want to interpolate the distortions
    ; INTERPOLATE uses the pixel coordinates of the input array to define the coordinate
    ; grid. So we need to define the coordinates based on the original pixels of the 
    ; original displacement array (which has now been extended)
    ; so first we find which pixel in the Rx/
    extend_size               = size(Rx)
    extend_start_x            = interpol(findgen(extend_size[1]),rx[*,0],grid_start_x)
    extend_start_y            = interpol(findgen(extend_size[2]),ry[0,*],grid_start_y)

    output_coord_x            = (findgen(scene_x) - grid_start_x)/grid_spacing_x + extend_start_x
    output_coord_y            = (findgen(scene_y) - grid_start_y)/grid_spacing_y + extend_start_y
        
    distortions_interpolate_x = interpolate(Px, output_coord_x, output_coord_y, cubic=-0.5, /grid)
    distortions_interpolate_y = interpolate(Py, output_coord_x, output_coord_y, cubic=-0.5, /grid)
    
    ans[*,*,0]                = distortions_interpolate_x
    ans[*,*,1]                = distortions_interpolate_y

endelse

return, ans
end

;
; **********************************************************
; ******************** FUNCTION: apod_mask  *******************
; **********************************************************
;+
; PURPOSE:
;  create an apodization mask
;
; INPUTS:
;  nx, ny: size of apodization window to create
;
; KEYWORD PARAMETERS:
;  taper_percent: what at what percentage of the array size (on each side)
;                    should the apodization roll off begin
; OUTPUTS:
;
; NOTES:
;   Optimal apodization percentage seems to be 0.08 for the Blackman window,
;     and 0.09 for the Hann window
; MODIFICATION HISTORY:
;  26 October 1990: Originally written by Phil Wiborg
;  August 2021: K. Reardon - cleaned up and modernized
;-
function	apod_mask, nx, ny, taper_percent=taper_percent

IF N_ELEMENTS(taper_percent) EQ 0 THEN taper_percent=0.08
IF (taper_percent LT 0) OR (taper_percent GT 1) THEN taper_percent=0.08

; makes a Hanning window with a width given by the edge percentage value
;edge_percentage = 10 ; percent
; mm=myhanning(nx,ny,10)
; but this eliminates the apodization filter, setting all values to unity.
; what? why?
;mm = fltarr(nx,ny) + 1
;mm(*)=1
; if we wanted no apodization, we can just use a value of zero for the taper_percent
;mm = apod_window(nx, ny, 0.0, 0.0, 2)

;window_type = 0   ; Hann window
window_type = 2   ; Blackman window

mm = apod_window(nx, ny, taper_percent, taper_percent, window_type )

return, mm
end
;
; **********************************************************
; ******************** FUNCTION: make_lowpass_filt  *******************
; **********************************************************
;+
; PURPOSE:
;  This function generates a low pass filter in FFT space
;
; INPUTS:
;  nx, ny = size of kernel subarray to be used
;
; KEYWORD PARAMETERS:
;  nxy_div = the normalization factor that determines how much the 
;              the high frequencies are attenuated
;              it appears that values below 6 cause further underestimation of the 
;              subfield displacements
;
; OUTPUTS:
;
; MODIFICATION HISTORY:
;  26 October 1990: Originally written by Phil Wiborg
;  August 2021: K. Reardon - cleaned up and modernized
;-
function	make_lowpass_filt, nx, ny, nxy_div=nxy_div

IF NOT keyword_set(nxy_div) THEN nxy_div = 6
divisor_limit = 10

x_idx = findgen (FIX(nx/2))
if (nx mod 2) EQ 1 then x_idx = [x_idx, x(nx/2-1), rotate(x_idx,2)] else x_idx = [x_idx, rotate(x_idx,2)]

x = exp(-(x_idx/(ROUND(nx/nxy_div) > divisor_limit))^2)

y_idx = findgen (FIX(ny/2))
if (ny mod 2) EQ 1 then y_idx = [y_idx, y_idx(ny/2-1), rotate(y_idx,2)] else y_idx = [y_idx, rotate(y_idx,2)]
y = exp(-(y_idx/(ROUND(ny/nxy_div) > divisor_limit))^2)

; matrix multiplication to make a FFT filer, with peak values in the corners.
mm = x # y

return, mm
end
;
; **********************************************************
; ******************** FUNCTION: doref  *******************
; **********************************************************
;+
; PURPOSE:
;  This function sets up reference windows
;
; INPUTS:
;  ref:   a 2-dimensional array (L x M) containing the reference 
;             against which the scene should be registered
;
; KEYWORD PARAMETERS:
;  subfield_images = the array of image subfields, as cutout from the full array
;  destr_info: 
; OUTPUTS:
;  subfield_fftconj = array of conjugate of FFT for all subfields
;
; MODIFICATION HISTORY:
;  August 2021: K. Reardon - cleaned up and modernized
;-
FUNCTION	doref, ref, apod_mask, subfield_images=subfield_images, destr_info=destr_info 	

; kx,ky = size of subfield box to calculate over
; cpx,cpy = number of subfield to calculate
; so this is an array of all subfields

subfield_fftconj = complexarr (destr_info.kx, destr_info.ky, destr_info.cpx, destr_info.cpy)
subfield_images  = fltarr (destr_info.kx, destr_info.ky, destr_info.cpx, destr_info.cpy)

; number of pixels in subfield
nelz   = destr_info.kx * destr_info.ky
xcoord = REBIN(FINDGEN(destr_info.kx,1), destr_info.kx, destr_info.ky)
ycoord = REBIN(FINDGEN(1,destr_info.ky), destr_info.kx, destr_info.ky)

sub_strt_y = destr_info.by
sub_end_y  = sub_strt_y + destr_info.ky - 1

for j = 0, destr_info.cpy-1 do begin

    sub_strt_x = destr_info.bx
    sub_end_x  = sub_strt_x + destr_info.kx - 1

    for i = 0, destr_info.cpx-1 do begin
        ; cut out subfield from reference

        ; instead of incrementing the subarray positions, we should just 
        ; take the reference positions from the predefined control points
        ; and calculate the subarray around those coordinates
        ; This will be more flexible going forward, especially considering 
        ; the possibility of irregular sampling
        sub_strt_x  = destr_info.rcps[0,i,j] - destr_info.kx/2
        sub_end_x   = sub_strt_x + destr_info.kx - 1

        sub_strt_y  = destr_info.rcps[1,i,j] - destr_info.ky/2
        sub_end_y   = sub_strt_y + destr_info.ky - 1

    	z = ref(sub_strt_x:sub_end_x, sub_strt_y:sub_end_y)

        CASE destr_info.subfield_correction OF
            'mean_subtraction': BEGIN
                IF (destr_info.debug GE 2) and N_ELEMENTS(scheme_printed) EQ 0 THEN BEGIN
                    print,'Subtracting mean value from subfields'
                    scheme_printed = 1
                ENDIF
                z -= MEAN(z)
              END

            'plane_subtraction': BEGIN
                IF (destr_info.debug GE 2) and N_ELEMENTS(scheme_printed) EQ 0 THEN BEGIN
                    print,'Subtracting fitted plane from subfields'
                    scheme_printed = 1
                ENDIF
                ; subtracts a fitted plane from the subfield (to remove intensity gradients)
                ; there are various methods to do this
                fit_degree = 1
	        z = z - sfit(z, fit_degree)
	        ;z=z-sfit(z, fit_degree,max_degree=1)
	        ;z_fit_param = planefit(xcoord, ycoord, z, 0.0, z_fit)
	        ;z = z - z_fit   
              END
            'surface_subtraction': BEGIN
                IF (destr_info.debug GE 2) and N_ELEMENTS(scheme_printed) EQ 0 THEN BEGIN
                    print,'Subtracting fitted 2nd order plynomial surface from subfields'
                    scheme_printed = 1
                ENDIF
                ; subtracts a fitted plane from the subfield (to remove intensity gradients)
                ; there are various methods to do this
                fit_degree = 2
	        z = z - sfit(z, fit_degree)
              END
            ELSE: z = z
        ENDCASE

        subfield_images(*,*,i,j) = z
	    
        ; stores the complex conjugate of the FFT of each reference subfield 
        ;    (for later calculation of the cross correlation with the target subfield)
        subfield_fftconj(*,*,i,j) = CONJ (FFT (z * apod_mask, -1))

        ; commands for zero-padding the array before the fft
        ;subfield_zeropad = fltarr(destr_info.kx*2, destr_info.ky*2)
        ;subfield_zeropad[destr_info.kx/2, destr_info.ky/2] = z * apod_mask
        ;subfield_fftconj(*,*,i,j) = CONJ (FFT (subfield_zeropad, -1))

        ; increments the starting point of the next subfield (in steps of kx)
	sub_strt_x = sub_strt_x + destr_info.spacing_x
	sub_end_x  = sub_strt_x + destr_info.kx - 1

    endfor ; loop over I

    sub_strt_y = sub_strt_y + destr_info.spacing_x
    sub_end_y  = sub_strt_y + destr_info.ky - 1

endfor

return, subfield_fftconj

end
;

; **********************************************************
; ******************** FUNCTION: cploc  *******************
; **********************************************************
;+
; PURPOSE:
;  This function sets up reference windows
;
; INPUTS:
;  scene:   a 2-dimensional array (L x M) containing the reference 
;             against which the scene should be registered
;
; KEYWORD PARAMETERS:
;  subfield_images = the array of image subfields, as cutout from the full array
;  destr_info: 
; OUTPUTS:
;  subfield_fftconj = array of conjugate of FFT for all subfields
;
; MODIFICATION HISTORY:
;  August 2021: K. Reardon - cleaned up and modernized
;-
FUNCTION	cploc, scene, subfield_fftconj, apod_mask, lowpass_filter, destr_info=destr_info	; locate control points

; output pixel values at each control point
subfield_offsets = fltarr(2, destr_info.cpx, destr_info.cpy)	; gets the results

;number of pixel elements in subfield
nels = destr_info.kx*destr_info.ky

xcoord = REBIN(FINDGEN(destr_info.kx,1), destr_info.kx, destr_info.ky)
ycoord = REBIN(FINDGEN(1,destr_info.ky), destr_info.kx, destr_info.ky)

sub_strt_y = destr_info.by
sub_end_y  = sub_strt_y + destr_info.ky - 1

subfield_crosscorrel_all = FLTARR(destr_info.kx, destr_info.ky, destr_info.cpx, destr_info.cpy)

for j = 0, destr_info.cpy-1 do begin

    sub_strt_x = destr_info.bx
    sub_end_x  = sub_strt_x + destr_info.kx - 1

    for i = 0, destr_info.cpx-1 do begin
        ; cut out subfield from reference

        ; instead of incrementing the subarray positions, we should just 
        ; take the reference positions from the predefined control points
        ; and calculate the subarray around those coordinates
        ; This will be more flexible going forward, especially considering 
        ; the possibility of irregular sampling
        sub_strt_x  = destr_info.rcps[0,i,j] - destr_info.kx/2
        sub_end_x   = sub_strt_x + destr_info.kx - 1

        sub_strt_y  = destr_info.rcps[1,i,j] - destr_info.ky/2
        sub_end_y   = sub_strt_y + destr_info.ky - 1

        scene_subfield = scene[sub_strt_x:sub_end_x, sub_strt_y:sub_end_y]

        CASE destr_info.subfield_correction OF
            'mean_subtraction': scene_subfield -= MEAN(scene_subfield)

            'plane_subtraction': BEGIN
                ; subtracts a fitted plane from the subfield (to remove intensity gradients)
                ; there are various methods to do this
                fit_degree = 1
	        scene_subfield = scene_subfield - sfit(scene_subfield, fit_degree)
	        ; scene_subfield = scene_subfield - sfit(scene_subfield, fit_degree,max_degree=1)
	        ;z_fit_param = planefit(xcoord, ycoord, scene_subfield, 0.0, subfield_plane)
	        ; scene_subfield = scene_subfield - subfield_plane   
              END
            'surface_subtraction': BEGIN
                ; subtracts a fitted plane from the subfield (to remove intensity gradients)
                ; there are various methods to do this
                fit_degree = 2
	        scene_subfield = scene_subfield - sfit(scene_subfield, fit_degree)
              END
            ELSE: scene_subfield = scene_subfield
        ENDCASE

        ; 1) take FFT of subfield
        ; 2) multiply it by complex conjugate of FFT of reference subfield
        ; 3) multiply the product by a low-pass filter      
        ; 4) take the absolute values (Fourier power)
        ; 5) shift 2-D power to center peak in array

        scene_subfield_fft      = fft(scene_subfield * apod_mask,-1)
        subfield_crosspower     = scene_subfield_fft * subfield_fftconj[*,*,i,j] * lowpass_filter

        ; zero pad the array before the fft
        ;subfield_zeropad = fltarr(destr_info.kx*2, destr_info.ky*2)
        ;subfield_zeropad[destr_info.kx/2, destr_info.ky/2] = scene_subfield * apod_mask
        ;scene_subfield_fft      = fft(subfield_zeropad,-1)
        ;subfield_crosspower     = scene_subfield_fft * subfield_fftconj[*,*,i,j] 
        ;subfield_crosspower     = subfield_crosspower[destr_info.kx/2:destr_info.kx*3/2-1, destr_info.ky/2:destr_info.ky*3/2-1]
        ;subfield_crosspower     = subfield_crosspower * lowpass_filter

        subfield_crosscorrel    = ABS(fft(subfield_crosspower, 1))
        subfield_crosscorrel    = SHIFT(subfield_crosscorrel, destr_info.kx/2,destr_info.ky/2)
	;cc = shift(abs(fft(fft(ss,-1)*w(*,*,i,j)*smou,1)),destr_info.wx/2,destr_info.wy/2)
        subfield_crosscorrel_all[*,*,i,j] = subfield_crosscorrel

       ; find peak of power and pixel location 
        mx = max (subfield_crosscorrel, loc)
    
    ;	simple maximum location - integer pixel accuracy
        ccsz = size (subfield_crosscorrel)
        xmax = loc mod ccsz(1)
        ymax = loc/ccsz(1)
    
        CASE destr_info.max_fit_method OF
        'niblack_simple' : BEGIN
            ;  a more complicated interpolation
            ; (from Niblack, W: An Introduction to Digital Image Processing, p 139.)
            IF (destr_info.debug GE 2) and N_ELEMENTS(scheme_printed) EQ 0 THEN BEGIN
                print,'Interpolating with 1st order Niblack scheme'
                scheme_printed = 1
            ENDIF

            ; test to see if identified maximum is away from the edges of the array
            if (xmax * ymax gt 0) and (xmax lt (ccsz(1)-1)) and (ymax lt (ccsz(2)-1)) then begin
                ; Sep 91 phw	try including more points in interpolations
                denom = mx*2 - subfield_crosscorrel[xmax-1,ymax] - subfield_crosscorrel[xmax+1,ymax]
                xfra = (xmax-0.5) + (mx - subfield_crosscorrel[xmax-1,ymax])/denom
    
                denom = mx*2 - subfield_crosscorrel[xmax,ymax-1] - subfield_crosscorrel[xmax,ymax+1]
                yfra = (ymax-0.5) + (mx - subfield_crosscorrel[xmax,ymax-1])/denom

                xmax=xfra
                ymax=yfra
            endif
          END

        'niblack_complex' : BEGIN
            ;A more complicated interpolation technique to find the maximum of an array.
            ;from, Niblack, W., "An Introduction to Digital Image Processing", p 139.
            IF (destr_info.debug GE 2) and N_ELEMENTS(scheme_printed) EQ 0 THEN BEGIN
                print,'Interpolating with 2nd order Niblack scheme'
                scheme_printed = 1
            ENDIF

            if (xmax*ymax gt 0) and (xmax lt (ccsz[1]-1)) and (ymax lt (ccsz[2]-1)) then begin
                a2 = (subfield_crosscorrel[xmax+1,ymax]    - subfield_crosscorrel[xmax-1,ymax])/2.
                a3 = (subfield_crosscorrel[xmax+1,ymax]/2. - subfield_crosscorrel[xmax,ymax]        + subfield_crosscorrel[xmax-1,ymax]/2.)
                a4 = (subfield_crosscorrel[xmax,ymax+1]    - subfield_crosscorrel[xmax,ymax-1])/2.
                a5 = (subfield_crosscorrel[xmax,ymax+1]/2. - subfield_crosscorrel[xmax,ymax]        + subfield_crosscorrel[xmax,ymax-1]/2.)
                a6 = (subfield_crosscorrel[xmax+1,ymax+1]  - subfield_crosscorrel[xmax+1,ymax-1] - $
                      subfield_crosscorrel[xmax-1,ymax+1]  + subfield_crosscorrel[xmax-1,ymax-1])/4.

                xdif = (2*a2*a5 - a4*a6)/(a6^2-4*a3*a5)
                ydif = (2*a3*a4 - a2*a6)/(a6^2-4*a3*a5)
                xmax = xmax + xdif
                ymax = ymax + ydif
            endif
          END
        'gaussian_fit' : BEGIN
            IF (destr_info.debug GE 2) and N_ELEMENTS(scheme_printed) EQ 0 THEN BEGIN
                print,'Interpolating with Gaussian fit'
                scheme_printed = 1
            ENDIF
            offs = 7
            if (xmax gt offs) and (ymax GT offs) and (xmax lt (ccsz(1)-offs-1)) and (ymax lt (ccsz(2)-offs-1)) then begin
                fitted_gauss = gauss2dfit(subfield_crosscorrel[xmax-offs:xmax+offs,ymax-offs:ymax+offs]^2,gauss_coeff)
                xmax = xmax - offs + gauss_coeff[4]
                ymax = ymax - offs + gauss_coeff[5]
            endif
          END
        'surface_fit' : BEGIN
            IF (destr_info.debug GE 2) and N_ELEMENTS(scheme_printed) EQ 0 THEN BEGIN
                print,'Interpolating with polynomial fit'
                scheme_printed = 1
            ENDIF
            boxsz = 4
            if (xmax gt boxsz) and (ymax GT boxsz) and (xmax lt (ccsz(1)-boxsz-1)) and (ymax lt (ccsz(2)-boxsz-1)) then begin
                fitted_surf = sfit(subfield_crosscorrel[xmax-boxsz:xmax+ boxsz,ymax-boxsz:ymax+ boxsz],2,kx=fitted_max_coeff)
                a0 = fitted_max_coeff[0]
                a1 = fitted_max_coeff[1]
                a2 = fitted_max_coeff[2]
                a3 = fitted_max_coeff[3]
                a4 = fitted_max_coeff[4]
                a5 = fitted_max_coeff[5]

                xdif = (2*a4*a1 - a5*a2)/(a5^2 - 4*a3*a4)
                ydif = (2*a3*a2 - a5*a1)/(a5^2 - 4*a3*a4)

                xmax = xmax + xdif
                ymax = ymax + ydif                 
            endif
          END
        ELSE: BEGIN
            xmax = xmax
            ymax = ymax
          END
        ENDCASE

        subfield_offsets(0,i,j) = sub_strt_x + xmax
        subfield_offsets(1,i,j) = sub_strt_y + ymax
    
        if destr_info.do_plots GE 1 then begin
;            plots, reform(subfield_offsets(*,i,j),2,1), /dev, psym=1,col=5e4

;            plots,(subfield_offsets(0,i,j) - destr_info.rcps[0,i,j])[0] * 5 + (destr_info.rcps[0,i,j]),$
;                  (subfield_offsets(1,i,j) - destr_info.rcps[1,i,j])[0] * 5 + (destr_info.rcps[1,i,j]),/dev,col=5e4,th=2,psym=1,syms=0.5
            
            vector_scale = MIN([destr_info.spacing_x,destr_info.spacing_y]) / 3.
            plots,(subfield_offsets(0,i,j) - destr_info.rcps[0,i,j])[0] * [0, vector_scale] + (destr_info.rcps[0,i,j]),$
                  (subfield_offsets(1,i,j) - destr_info.rcps[1,i,j])[0] * [0, vector_scale] + (destr_info.rcps[1,i,j]),/dev,col=200,th=2
            wait,0.	; flushes plots
        endif
    
        ; increments the starting point of the next subfield (in steps of kx)
        sub_strt_x = sub_strt_x + destr_info.spacing_x
        sub_end_x  = sub_strt_x + destr_info.kx - 1
    endfor

    sub_strt_y = sub_strt_y + destr_info.spacing_y
    sub_end_y  = sub_strt_y + destr_info.ky - 1
endfor

return, subfield_offsets

end
;	

; **********************************************************
; ********************* FUNCTION: doreg  *******************
; **********************************************************
;+
; PURPOSE:
;  This function sets up reference windows
;
; INPUTS:
; scene: fltarr(nx,ny)     -  image to be destretched
; rdisp: fltarr(2,cpx,cpy) -  reference control points
; disp:  fltarr(2,cpx,cpy) -  actual displaced position of control points
;
; KEYWORD PARAMETERS:
;  subfield_images = the array of image subfields, as cutout from the full array
;  destr_info: 
; OUTPUTS:
;  destretched_scene = array after application of displacements
;
; MODIFICATION HISTORY:
;  August 2021: K. Reardon - cleaned up and modernized
;-

function doreg, scene, rdisp, disp	; destretch by re-sample of scene

;	old method
;degree = 1			; linear (1st degree) interpolation works best
;cx = fltarr (degree +1, degree +1,/nozero)
;cy = fltarr (degree +1, degree +1,/nozero)
;polywarp, disp(0,*,*), disp(1,*,*), rdisp(0,*,*), rdisp(1,*,*), degree, cx, cy
; destretched_scene = poly_2d (scene, cx, cy, 1)

; these are newer methods for doing the interpolation

; this is a procedure in this file that interpolates the measured shifts at a subgrid of 
;     control points onto a full array equal to the input scene size (which comes from destr_info)
;     original_version=1 uses the original B-spline algorithm
;     original_version=0 uses IDL's built-in INTERPOLATE command
xy = rescale_distorion_map (rdisp, disp, destr_info=destr_info, original_version=0)

; this is a procedure in this file, but there is also a bilin.pro procedure from Ray Sterner
; destretched_scene = bilin (scene, xy, bilinear=0)	; nearest neighbor interpolation
destretched_scene = bilin (scene, xy, bilinear=1)	; bi-linear interpolation

; this is a builtin procedure with IDL ($IDL_RID/lib/bilinear.pro)
; http://www.harrisgeospatial.com/docs/bilinear.html
; destretched_scene = bilinear (scene, xy(*,*,0), xy(*,*,1))	; from userLib: may not work!

; What about using kgrid2d? 
; http://www.harrisgeospatial.com/docs/krig2d.html#

return, destretched_scene
end
;

; **********************************************************
; ********************  PRO: setup  *******************
; **********************************************************;+
; PURPOSE:
;  This function performs some initial checks and evaluation of input arrays
;
; INPUTS:
;  scene: a 2 or 3-dimensional array (L x M x N) containing the images to be destretched
;  ref:   a 2-dimensional array (L x M) containing the reference 
;             against which the scene should be registered
;  kernel: a 2-dimensional array (S x T) defining the size of the cross-correlation subfields
;
; KEYWORD PARAMETERS:
;  destr_info: shared variable containing destretching parameters and control points
;
; OUTPUTS:
;
; MODIFICATION HISTORY:
;  26 October 1990: Originally written by Phil Wiborg
;  August 2021: K. Reardon - cleaned up and modernized
;-
PRO	setup, scene, ref, kernel, destr_info=destr_info	; initialize cps & reg

IF N_ELEMENTS(destr_info) EQ 0 THEN BEGIN
    destr_info = define_destr_info()
    destr_info.debug = 3
ENDIF

; instead of yet another procedure, just move the checks from val into this procedure
;val, scene, ref, kernel

scene_size = size(scene)
reference_size = size(ref)
kernel_size = size(kernel)
error_flag = 0

if (scene_size(0) ne 2) and (scene_size(0) ne 3) then begin
    print, "ERROR: Input scene(s) must be in a 2-D or 3-D array"
    error_flag = error_flag + 1
endif

if reference_size(0) ne 2 then begin
    print, "ERROR: Destretching reference must be a 2-D array"
    error_flag = error_flag + 1
    endif

if (scene_size(1) ne reference_size(1)) or (scene_size(2) ne reference_size(2)) then begin
    print, "ERROR: Both the x and y dimensions of the input scene must match those of the reference"
    error_flag = error_flag + 1
    endif

if kernel_size(0) ne 2 then begin
    print, "ERROR: Destretching kernel must be 2-D array"
    error_flag = error_flag + 1
    endif

if error_flag GT 0 then stop, "ERROR: Quitting - too many errors"

; determine dimensions of input scene - is it 2-D or 3-D?
; if it is only a 2-D array add a dummy third dimension for consistency
IF scene_size[0] eq 2 THEN BEGIN
    scene = reform (scene, scene_size[1], scene_size[2], 1)
    scene_size = size(scene)
ENDIF

destr_info.kx         = kernel_size[1]
destr_info.ky         = kernel_size[2]

destr_info.ref_sz_x   = reference_size[1]
destr_info.ref_sz_y   = reference_size[2]

destr_info.scene_sz_x = scene_size[1]
destr_info.scene_sz_y = scene_size[2]
destr_info.scene_sz_z = scene_size[3]


return

end
; **********************************************************
; ********************  END: setup  *******************
; **********************************************************;+



; **********************************************************
; ********************  FUNCTION: repair  ******************
; **********************************************************
; PURPOSE:
;  This function finds & fixes 'bad' control points, where
;      a bad point is one where the calculated displacement is
;      too far from the reference position
;
; INPUTS:
;  refrnc_cntl_pts:   fltarr(2, nx, ny) an array containing the reference coordinates
;  disp:  fltarr(2, nx, ny, nz) an array of displacements to be checked
;
; KEYWORD PARAMETERS:
;  max_displacement: the allowed displacement limit, given in terms of a fraction of
;                      the kernel size
;  destr_info: shared variable containing destretching parameters and control points
;
; OUTPUTS:
;
; MODIFICATION HISTORY:
;  26 October 1990: Originally written by Phil Wiborg
;  August 2021: K. Reardon - cleaned up and modernized
;-

FUNCTION	repair, refrnc_cntl_pts, disp, max_displacement=max_displacement, destr_info=destr_info

IF N_ELEMENTS(max_displacement) EQ 0 THEN max_displacement=0.25

IF (max_displacement LE 0) OR (max_displacement GE 1) THEN max_displacement=0.25

; explicitly determine size of displacement array, including determining whether 
;     it is has a fourth dimension
sz = size (disp)
nd = sz[1]
IF nd NE 2 THEN print,'displacement array should be an array of the form [2,nx,ny,nf]'
nx = sz[2]
ny = sz[3]
if (sz[0] eq 4) then nf = sz[4] else nf = 1

; the maximum allowable offset from the control point is the larger of the
; kernel size in the x and y directions
; this is done for simplicity, considering that kernels are often essentially square,
; but it would be straightforward to perform separate checks for both kernel axes
limit = (max ([destr_info.kx, destr_info.ky]) * max_displacement)

; make a copy of the disp array into which to write any corrected values
disp_repaired = disp
for frm=0, nf-1 do begin
    ; select out one set of determined offsets
    disp_1frame = reform(disp[*,*,*,frm])
    offsets_calc = disp_1frame - refrnc_cntl_pts
    offsets_calc_mag = offsets_calc[0,*,*]^2 + offsets_calc[1,*,*]^2

    ; list of bad coordinates in this frame
    ; searches for where the quadratic sum of the offsets or the individual offsets 
    ; in either of the two x and y axes is greater than the provided limit
    bad_offset = where ((offsets_calc_mag GT limit^2) OR $
                        (offsets_calc[0,*,*] GT limit) OR (offsets_calc[1,*,*] GT limit), count)
    
    offsets_median_x = MEDIAN(REFORM(offsets_calc[0,*,*]),3)
    offsets_median_y = MEDIAN(REFORM(offsets_calc[1,*,*]),3)

    for i = 0, count-1 do begin
        x = bad_offset(i) mod nx
        y = bad_offset(i)/nx
        if destr_info.do_plots GE 1 then begin
            symsize_scaled = (SQRT(MIN([destr_info.spacing_x,destr_info.spacing_y]) / 24)) > 0.5 < 2
            plots, reform (disp_1frame(*,x,y),2), /dev, psym=4, col=5e4, thick=1,syms=symsize_scaled
        endif

;        Repairs disp_repaired @ (x,y). 
;        The best thing would be to replace bad points with some local
;        mean or median. But there are many cases to consider:
;        is the point needing repair at one of 4 corners, on
;        borders, or in interior?  Also, how many other neighbors
;        are also bad points?  
;        Therefore, we take the easy way
;        out (for now) and say the 'repair' is 'no warping'. 
;	disp_repaired (*,x,y,frm) = refrnc_cntl_pts(*,x,y)

; if we just run a median filter over the whole array, we can replace the bad displacement in
;     one position with the median of the surrounding points. This might have some issues at the
;     edges but might otherwise be a small improvement of setting the offsets to zero
	disp_repaired(0,x,y,frm) = refrnc_cntl_pts(0,x,y) + offsets_median_x[x,y]
	disp_repaired(1,x,y,frm) = refrnc_cntl_pts(1,x,y) + offsets_median_y[x,y]
        
    endfor

endfor

return, disp_repaired
end
;
; ****************************************************************
; ********************* FUNCTION: nee reg    *********************
; ********************* FUNCTION: destretch  *********************
; ****************************************************************
function destretch, scene, ref, kernel, apply_scene=apply_scene,$
              rdisp=rdisp, disp=disp,$
              crosscor_fit=crosscor_fit, lowpass_factor=lowpass_factor , apod_percent=apod_percent, $
              border_offset=border_offset, spacing_ratio=spacing_ratio,$
              debug_level=debug_level, plot_level=plot_level
; register scene(s) with respect to ref, using kernel size

; scene is (nx,xy) or (nx,ny,nf)	scene(s) to be registered
; ref is (nx,ny)			reference scene
; kernel is (kx,ky)			conveys the size of the kernel
;					and is otherwise unused

; returns destretched scene: (nx,ny[,nf])

IF N_ELEMENTS(border_offset)  EQ 0 THEN border_offset = 2
IF N_ELEMENTS(spacing_ratio)  EQ 0 THEN spacing_ratio = 0.5
IF NOT KEYWORD_SET(apply_scene)    THEN apply_scene=scene
IF NOT KEYWORD_SET(crosscor_fit)   THEN crosscor_fit ='niblack_complex'
IF NOT KEYWORD_SET(lowpass_factor) THEN lowpass_factor = 6
;   Optimal apodization percentage seems to be 0.08 for the Blackman window
IF N_ELEMENTS(apod_percent)   EQ 0 THEN apod_percent = 0.08

debug = 0
; pretty basic tests to make sure arrays are of appropriate dimensions and matching sizes
; this program also sets the values in destr_info for the sizes of the kernel, scene, and reference
;     destr_info.[kx,ky], destr_info.[ref_sz_x,ref_sz_y], destr_info.[scene_sz_x,scene_sz_y,scene_sz_z]
setup, scene, ref, kernel, destr_info=destr_info

; set level of debugging output
IF N_ELEMENTS(debug_level) EQ 0 THEN destr_info.debug = 3    ELSE destr_info.debug = debug_level
; set amount of output plotting
IF N_ELEMENTS(plot_level)  EQ 0 THEN destr_info.do_plots = 1 ELSE destr_info.do_plots = plot_level

; determine border and control point matrix size
rdisp = destr_control_points(ref, kernel, destr_info=destr_info, border_offset=border_offset, spacing_ratio=spacing_ratio)

; make an apodization mask
; use taper_percent=0.0 to generate a mask of all 1's = no apodization
apod_subarray = apod_mask (destr_info.kx, destr_info.ky, taper_percent=apod_percent)

; make an FFT frequency (low pass) filter (high in the corners)
;print,wx,wy,kx,ky,bx,by,cpx,cpy
lowpass_subarray = make_lowpass_filt(destr_info.kx * 2, destr_info.ky * 2, nxy_div=lowpass_factor)

; select different options for correcting subfields - leave only one option uncommented.
;   -- the mean subtraction is important before performing the apodization and FFT
;   -- the plane subtraction removes the mean but also trends, which can be useful
; in regions with gradients, like near the limb or sunspots, for example
;   -- should this option be user selectable at some point? maybe...
;destr_info.subfield_correction = 'mean_subtraction'
;destr_info.subfield_correction = 'plane_subtraction'
destr_info.subfield_correction = 'surface_subtraction'
;destr_info.subfield_correction = 'none'

destr_info.max_fit_method = crosscor_fit
;destr_info.max_fit_method = 'niblack_simple'
;destr_info.max_fit_method = 'niblack_complex'
;destr_info.max_fit_method = 'gaussian_fit'

; create an array of subfields from the reference image
; these are actually the complex conjugate of the FFT of each subfield
subfield_fftconj = doref(ref, apod_subarray, subfield_images=subfield_images, destr_info=destr_info)

; define the size of the input and output target arrays (scene)
destretched_scene = fltarr (destr_info.scene_sz_x, destr_info.scene_sz_y, destr_info.scene_sz_z,/nozero)

; compute control point locations

; step through each frame
for frm = 0, destr_info.scene_sz_z-1 do begin
    if destr_info.do_plots GE 1 then begin
        ; this plots a box showing the size of the destretch kernel centered on the lower left control point
        plots, rdisp[0,0,0] + destr_info.kx * [0.5,0.5],   rdisp[1,0,0] + destr_info.ky*[-0.5,0.5],/dev,col=5e4,th=2
        plots, rdisp[0,0,0] + destr_info.kx * [-0.5,-0.5], rdisp[1,0,0] + destr_info.ky*[-0.5,0.5],/dev,col=5e4,th=2
        plots, rdisp[0,0,0] + destr_info.kx * [-0.5,0.5],  rdisp[1,0,0] + destr_info.ky*[0.5,0.5],/dev,col=5e4,th=2
        plots, rdisp[0,0,0] + destr_info.kx * [-0.5,0.5],  rdisp[1,0,0] + destr_info.ky*[-0.5,-0.5],/dev,col=5e4,th=2

        if destr_info.do_plots GE 2 then begin
            ; this plots a square symbol at the location of each control point
            symsize_scaled = (SQRT(MIN([destr_info.spacing_x,destr_info.spacing_y]) / 24)) > 0.5 < 2
	    plots, reform(rdisp,2,n_elements(rdisp)/2), /dev, psym=6, thick=1, symsize=symsize_scaled
        endif
	;wait,0.		; enables flushing of plot!
    endif

    ; one image from the scene, the array of reference subfields, the apodization mask, and the FFT low-pass filter
    disp = cploc (scene[*,*,frm], subfield_fftconj, apod_subarray, lowpass_subarray, destr_info=destr_info)

    disp = repair (rdisp, disp, destr_info=destr_info) ; optional repair
    IF destr_info.debug GE 2 THEN BEGIN
        rms = sqrt(total((rdisp - disp)^2)/n_elements(rdisp[0,*,*]))
        print, 'Displacement rms =', rms
    ENDIF

    x = doreg (apply_scene(*,*,frm), rdisp, disp)
    destretched_scene(*,*,frm) = x
;    subfield_fftconj = doref (x, apod_subarray); optional update of window
    endfor

undo	; undo setup

if destr_info.scene_sz_z eq 1 then begin
    scene             = reform (scene, destr_info.scene_sz_x, destr_info.scene_sz_y)
    destretched_scene = reform (destretched_scene, destr_info.scene_sz_x, destr_info.scene_sz_y)
    endif

return, destretched_scene
end
;
;	initialization routine

common reg_save, x_sys, y_sys

x_sys = !x & y_sys = !y

end
