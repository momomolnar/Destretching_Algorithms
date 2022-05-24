; *************************************************************************
; ********************  FUNCTION: destr_control_points  *******************
; ********************  FUNCTION: nee mkcps             *******************
; *************************************************************************

; this function defines a regularly spaced grid on control points, which are
;     the central pixel positions of each subfield for the destretching local 
;     offset determination.

function	destr_control_points, destr_info=destr_info, $
    border_offset=border_offset, spacing_ratio=spacing_ratio	; choose control point locations

; reference is passed simply to determine the size of the reference
; kernel gives kernel size

; All control info is passed between support programs via common block
;common reg_com, kx, ky,		$; kernel x,y size (16,16)
;		wx, wy,		$; border offset (32,32)
;		bx, by,		$; boundary x,y size (4,4)
;		cpx, cpy,	$; control point x,y size (7,4)
;		debug		; = 0	silent operation
				; = 1	TV_FLG
				; = 2	GRD_FLG
				; = 4	PRT_FLG

do_it_the_old_way = 0

IF N_ELEMENTS(border_offset) EQ 0 THEN border_offset = 0
IF N_ELEMENTS(spacing_ratio)  EQ 0 THEN spacing_ratio = 1

IF N_ELEMENTS(destr_info) EQ 0 THEN BEGIN
    destr_info = define_destr_info()
    destr_info.debug = 3
ENDIF

; determine the number of pixels in the kernel
;kernel_size   = (size(kernel))[1:2]
;destr_info.kx = kernel_size[0]
;destr_info.ky = kernel_size[1]
IF destr_info.debug GE 2 then print,'Kernel Size = ',STRTRIM(destr_info.kx,2), ' x ', STRTRIM(destr_info.ky,2)

; determine the number of pixels in the reference image
;reference_size = size(reference)
; the assumption is that the reference is a 2D array, so we only need the x- and y-dimensions
;reference_size = reference_size[1:2]
;destr_info.ref_sz_x  = reference_size[0]
;destr_info.ref_sz_y  = reference_size[1]
IF destr_info.debug GE 2 then print,'Scene Size = ',STRTRIM(destr_info.ref_sz_x,2), ' x ', STRTRIM(destr_info.ref_sz_y,2)

; [wx,wy] define the size of a border around the edge of the image, to add an additional 
;     buffer area in which to avoid placing the control points.
; The border_offset input variable defines this border area in relation to the kernel size,
;     but maybe it's better to define it as an absolute number of pixels?
destr_info.border_x    = fix(border_offset)
destr_info.border_y    = fix(border_offset)
; make sure [border_x,border_y] is divisible by 2
if destr_info.border_x mod 2 eq 1 then destr_info.border_x = destr_info.border_x + 1
if destr_info.border_y mod 2 eq 1 then destr_info.border_y = destr_info.border_y + 1
IF destr_info.debug GE 2 then print,'Border Size = ',STRTRIM(destr_info.border_x,2), ' x ', STRTRIM(destr_info.border_y,2)

if keyword_set(do_it_the_old_way) then begin

    ; Define how many subfields will fit across image
    ; (-wx + kx) = -(wx - kx) = -(kx * border_offset - kx) = -(border_offset - 1)*kx
    ; it seems like (reference_size[nn] - wx)/kx would be sufficient
    destr_info.cpx = (destr_info.ref_sz_x - destr_info.border_x + destr_info.kx)/destr_info.kx
    destr_info.cpy = (destr_info.ref_sz_y - destr_info.border_y + destr_info.ky)/destr_info.ky
    IF destr_info.debug GE 2 then print,'Number of Control Points = ',STRTRIM(destr_info.cpx,2), ' x ', STRTRIM(destr_info.cpy,2)

    ; Now define the starting pixels for the control points, taking into account
    ;     both the border and the width of the kernel
    ; allocate half of "leftover" pixels on either side of grid-point array 
    destr_info.bx = ((destr_info.ref_sz_x - destr_info.border_x + destr_info.kx) mod destr_info.kx)/2
    destr_info.by = ((destr_info.ref_sz_y - destr_info.border_y + destr_info.ky) mod destr_info.ky)/2
    IF destr_info.debug GE 2 then print,'Number of Border Pixels = ',STRTRIM(destr_info.bx,2), ' x ', STRTRIM(destr_info.by,2)
 
    ; setup pixel positions of grid points
    rcps = fltarr (2,cpx,cpy)
    ; ly = offset position (edge boundary) in y-direction
    ly = destr_info.by
    ; hy = starting position of grid - one full border offset in from boundary
    hy = ly + destr_info.border_y
    for j = 0, destr_info.cpy-1 do begin
        lx = destr_info.bx
        hx = lx + destr_info.border_x
        for i = 0, destr_info.cpx-1 do begin
            ; confusing, but its really a grid with the step size of the kernel, offset by half the width of the 
            ;   border limits, plus the boundary buffer
            ; i = 0 : rcps(0,i,j) = (lx + hx)/2. = (2 * bx + wx) / 2. = (bx + wx/2.)
            ; i = 1 : rcps(0,i,j) = ((bx + kx) + (bx + wx) + kx ) / 2. = (bx + wx/2. + kx)
            ; i = 2 : rcps(0,i,j) = ((bx + kx + kx) + bx + wx + kx + kx) / 2. = (bx + wx/2. + kx * 2)
            ; i = 3 : rcps(0,i,j) = ((bx + kx*3) + bx + wx + kx*3) / 2. = (bx + wx/2. + kx * 3)
            ; ...
            ; i = n : rcps(0,i,j) = ((bx + kx*n) + bx + wx + kx*n) / 2. = (bx + wx/2. + kx * n)
            rcps(0,i,j) = (lx + hx)/2
            rcps(1,i,j) = (ly + hy)/2
            lx = lx + kx
            hx = hx + kx
        endfor
        ly = ly + ky
        hy = hy + ky
    endfor

endif else begin

    ; the control points must start and end at least 1/2 kernel width away from the edges of the array
    ; So that means that the allowable range of pixels available for control points 
    ;     is reduced by (at minimum) one kernel width
    allowable_range_x = destr_info.ref_sz_x - destr_info.kx - (destr_info.border_x * 2)
    allowable_range_y = destr_info.ref_sz_y - destr_info.ky - (destr_info.border_y * 2)

    ; how far apart should the sub-array control points be placed, in units of the kernel width
    ; set the spacing between subarrays, making sure it is divisible by 2 (just because...)
    destr_info.spacing_x = fix(destr_info.kx * spacing_ratio)
    destr_info.spacing_y = fix(destr_info.ky * spacing_ratio)
    destr_info.spacing_x += destr_info.spacing_x mod 2
    destr_info.spacing_y += destr_info.spacing_y mod 2
    ;grid_spacing_x    = ROUND(destr_info.kx / spacing_ratio)
    ;grid_spacing_y    = ROUND (destr_info.ky / spacing_ratio)

    ; divide the number of allowable pixels by the control points, round down to nearest integer
    num_grid_x        = FIX(allowable_range_x / destr_info.spacing_x) + 1
    num_grid_y        = FIX(allowable_range_y / destr_info.spacing_y) + 1

    ; how far apart will the first and last control points be, in each axis
    total_range_x     = destr_info.spacing_x * (num_grid_x - 1)
    total_range_y     = destr_info.spacing_y * (num_grid_y - 1)
    ; the total range will be less than the maximum possible range, in most cases
    ; so allocate some of those extra pixels to each border
    destr_info.bx     = ROUND((allowable_range_x - total_range_x + destr_info.kx)/2.)
    destr_info.by     = ROUND((allowable_range_y - total_range_y + destr_info.ky)/2.)
    
    IF destr_info.debug GE 2 then print,'Number of Control Points = ',STRTRIM(num_grid_x,2), ' x ', STRTRIM(num_grid_y,2)
    IF destr_info.debug GE 2 then print,'Number of Border Pixels = ',STRTRIM(destr_info.bx,2), ' x ', STRTRIM(destr_info.by,2)
    IF destr_info.debug GE 3 THEN PRINT,'allowable range,grid spacing x, num grid x, total range x, start pos x', $
	                                   allowable_range_x,destr_info.spacing_x,num_grid_x,total_range_x,destr_info.bx

    rcps              = FLTARR(2, num_grid_x, num_grid_y)
    rcps[0,*,*]       = rebin(indgen(num_grid_x, 1) * destr_info.spacing_x + destr_info.bx, num_grid_x, num_grid_y) 
    rcps[1,*,*]       = rebin(indgen(1, num_grid_y) * destr_info.spacing_y + destr_info.by, num_grid_x, num_grid_y) 
    destr_info.cpx    = num_grid_x
    destr_info.cpy    = num_grid_y

endelse

; the size of the rcps array might have changed compared to the input destr_info, 
;     so we need to create a new structure, not just update the values in the old one

rcps_size = size(rcps)
destr_info_output = define_destr_info(rcps_size=rcps_size[1:3])

destr_info_output.kx           = destr_info.kx
destr_info_output.ky           = destr_info.ky
destr_info_output.spacing_x    = destr_info.spacing_x
destr_info_output.spacing_y    = destr_info.spacing_y
destr_info_output.border_x     = destr_info.border_x
destr_info_output.border_y     = destr_info.border_y
destr_info_output.wx           = destr_info.wx
destr_info_output.wy           = destr_info.wy
destr_info_output.bx           = destr_info.bx
destr_info_output.by           = destr_info.by
destr_info_output.cpx          = destr_info.cpx
destr_info_output.cpy          = destr_info.cpy
destr_info_output.ref_sz_x     = destr_info.ref_sz_x
destr_info_output.ref_sz_y     = destr_info.ref_sz_y
destr_info_output.scene_sz_x   = destr_info.scene_sz_x
destr_info_output.scene_sz_y   = destr_info.scene_sz_y
destr_info_output.scene_sz_z   = destr_info.scene_sz_z
destr_info_output.rcps         = rcps
destr_info_output.do_plots     = destr_info.do_plots
destr_info_output.debug        = destr_info.debug

destr_info = destr_info_output

IF destr_info.debug GE 2 THEN PRINT, 'Control points defined'

return, rcps

end
