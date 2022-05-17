
rot_angle    = 1
noise_level  = 0.05
cent         = 500
spacrat      = 0.5
ksz          = 64

add_noise    = 1

data_norm = data[*,*,1]/MEAN(data[*,*,1])
data_norm_rot = ROT(data_norm, rot_angle,1.0,cent,cent,/pivot,cubic=-0.5)

IF add_noise EQ 1 THEN BEGIN
    seed1 = systime(/sec)                            
    seed2 = systime(/sec)                            
    seed1 = (seed1 * 1e5 - LONG64(seed1 * 1e5)) * 1e6
    seed2 = (seed2 * 1e5 - LONG64(seed2 * 1e5)) * 1e6

    data_noise1 = (data_norm)     + (data_norm * randomn(seed1,[1000,1000]) * noise_level)
    data_noise2 = (data_norm_rot) + (data_norm_rot * randomn(seed2,[1000,1000]) * noise_level)
ENDIF ELSE BEGIN
    data_noise1 = data_norm
    data_noise2 = data_norm_rot
ENDELSE

erase
data_noise2_reg = destretch(data_noise2^2, data_noise1^2, bytarr(ksz, ksz),rdisp=rdisp,disp=disp,$
                spacing_ratio=spacrat,crosscor_fit='niblack_complex', lowpass_factor=6, apod_percent=0.08)


beta               = atan(rdisp[0,*,*] - cent, rdisp[1,*,*] - cent)                                                                                  
alpha              = rot_angle * !DTOR
distance           = sqrt(TOTAL((rdisp - cent)^2,1))
rot_shifts_calc_x  = (sin(beta + alpha) * distance) - (sin(beta) * distance)
rot_shifts_calc_y  = (cos(beta + alpha) * distance) - (cos(beta) * distance)
rot_shifts_destr_x = disp[0,*,*] - rdisp[0,*,*]
rot_shifts_destr_y = disp[1,*,*] - rdisp[1,*,*]            

print,median(rot_shifts_destr_x / rot_shifts_calc_x) , $
      median(rot_shifts_destr_y / rot_shifts_calc_y) 

end