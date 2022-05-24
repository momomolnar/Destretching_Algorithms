#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 16:09:55 2020

Destretching routines for removing optical defects following the
reg.pro by Phil Wiborg and Thomas Rimmele in reg.pro

@author: molnarad
"""

from asyncio import selector_events
import numpy as np
import scipy as sp
from scipy import signal as signal
from scipy.ndimage.interpolation import shift
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as pl
from time import time

class Destretch_params():
    """
    Class containing all the information about then

    """
    def __init__(self, kx, ky, wx, wy, bx, by, cpx, cpy, mf=0.08, ref_sz_x, ref_sz_y, 
                       scene_sz_x, scene_sz_y, subfield_correction, 
                       max_fit_method, do_plots, debug):
        self.kx = kx    # kernel size x
        self.ky = ky    # kernel size y
        self.wx = wx    # border offset in x direction
        self.wy = wy    # border offset in y direction
        self.bx = bx    # boundary size in x direction
        self.by = by    # boundary size in y direction
        self.cpx = cpx  # number of control points in x direction
        self.cpy = cpy  # number of control points in y direction
        self.mf = mf    # 
        self.ref_sz_x = ref_sz_x 
        self.ref_sz_y = ref_sz_y 
        self.scene_sz_x = scene_sz_x 
        self.scene_sz_y = scene_sz_y 
        self.subfield_correction = subfield_correction
        self.max_fit_method = max_fit_method
        self.do_plots = do_plots
        self.debug = debug

    def print_props(self):
        print("[kx, ky, wx, wy, bx, by, cpx, cpy, mf] are:")
        print(self.kx, self.ky, self.wx, self.wy, self.bx, self.by,
              self.cpx, self.cpy, self.mf)

def plot_cps(ax_object, d_info):
    """
    Plot the control points for the destretching on the destretched
    image

    Parameters
    ----------
    ax_object : matplotlib ax object
        Axis object to have the control points plotted on;
    d_info : class Destretch_params
        Destretch Parameters;

    Returns
    -------
    """
    return 0

def bilin_values_scene(scene, coords_new, destr_info, nearest_neighbor = False):
    """
    Bilinear interpolation (resampling)
    of the scene s at coordinates xy

    Parameters
    ----------
    scene : ndarray (nx, ny)
        Scene
    coords_new : ndarray (2, nx, ny)
        coordinates of the pixels of the output image
        on the input image (at which to interpolate the scene)
    destr_info: class Destretch_Params
        Destretch parameters

    Returns
    -------
    ans: ndarray (nx, ny)
        Bilinear interpolated (resampled) image at the xy locations

    """

    if nearest_neighbor == True:
        x = np.array(coords_new[:, :, 0] + .5, order="F")
        y = np.array(coords_new[:, :, 1] + .5, order="F")
        
        scene_interp = scene[x, y]

    else:
        x = np.array(coords_new[0, :, :], order="F")
        y = np.array(coords_new[1, :, :], order="F")

        # need to limit output coordinates so that interpolation calculations
        # don't go out of bounds (i.e. we add 1 to x and y coordinates below)
        x = np.clip(x, 0, x.shape[0]-2)
        y = np.clip(y, 0, y.shape[1]-2)

        x0 = x.astype(int)
        x1 = (x+1).astype(int)
        y0 = (y).astype(int)
        y1 = (y+1).astype(int)

        fx = x % 1.
        fy = y % 1.

        scene  = np.array(selector_events, order="F")
        scene_float = s.astype(scene)

        ss00 = scene_float[x0, y0]
        ss01 = scene_float[x0, y1]
        ssfx00 =                (scene_float[x1, y0] - ss00) * fx
        ssfy01 = (ss01 - ss00 + (scene_float[x1, y1] - ss01) * fx - ssfx00) * fy
        scene_interp  = ss00 + ssfx00 + ssfy01

    return scene_interp

def bilin_control_points(scene, rdisp, disp,
                         test=False):
    """
    Compute the coordinates of the pixels in the output images to be
    sampled from the input image (using Scipy.interpolate.RectBivariate).
    Interpolate the control point displacements to infer the sampling
    coordinates.

    Parameters
    ----------
    scene : ndarray (nx, ny)
        Image input
    rdisp : ndarray (kx, ky, 2)
        Reference coordinates of the control points.
    disp : ndarray (kx, ky, 2)
        Actual coordinates of the control points.

    Returns
    -------
    xy_grid : ndarray (2, nx, ny)
        Coordinates of the input image to be sampled for the output image

    """

    scene_nx = scene.shape[0]
    scene_ny = scene.shape[1]

    #compute the control points locations
    cp_x_coords = rdisp[0, :, 0]
    cp_y_coords = rdisp[1, 0, :]

    #compute the displacements

    xy_ref_coordinates1 = np.zeros((2, scene_nx, scene_ny), order="F")
    xy_ref_coordinates = np.zeros((2, scene_nx, scene_ny), order="F")




    xy_ref_coordinates[0, :, :] = [np.linspace(0, (scene_nx-1),
                                               num=scene_ny, dtype="int")
                                   for el in range(scene_nx)]
    xy_ref_coordinates[1, :, :] = [np.zeros(scene_ny, dtype="int")+el
                                   for el in range(scene_nx)]

    xy_ref_coordinates = np.swapaxes(xy_ref_coordinates, 1, 2)




    dd = disp - rdisp


    interp_x = RectBivariateSpline(cp_x_coords, cp_y_coords, dd[0, :, :])
    interp_y = RectBivariateSpline(cp_x_coords, cp_y_coords, dd[1, :, :])

    xy_grid = np.zeros((2, scene_nx, scene_ny))

    x_coords_output = np.linspace(0, scene_nx-1, num=scene_nx)
    y_coords_output = np.linspace(0, scene_ny-1, num=scene_ny)

    xy_grid[1, :, :] = 1. * interp_x.__call__(x_coords_output, y_coords_output,
                                         grid=True)
    xy_grid[0, :, :] = 1. *interp_y.__call__(x_coords_output, y_coords_output,
                                         grid=True)
    if test == True:
        im1 = pl.imshow(xy_grid[0, :, :])
        pl.colorbar(im1)
        pl.show()

        im2 = pl.imshow(xy_grid[1, :, :])
        pl.colorbar(im2)
        pl.show()

    xy_grid += xy_ref_coordinates

    return (xy_grid)

def bspline(scene, r, dd, d_info):
    """
    Destretch the scene using a B-spline

    Parameters
    ----------
    scene : TYPE
        Image to be destretched.
    r : TYPE
        reference displacements of control points
    dd : TYPE
        actual displacements of control points

    d_info: Destretch_class
        info about the destretching

    Returns
    -------
    ans : TYPE
        Destretched image

    """

    always = 1      # exterior control points drift with interior (best)
    #always = 0     ; exterior control points fixed by ref. displacements

    # a kludgery: increases magnitude of error, since
    # curve doesn't generally pass through the tie pts.
    d = (dd-r)*1.1 + r

    ds = r[0, 1, 0]-r[0, 0, 0]
    dt = r[1, 0, 1]-r[1, 0, 0]

    dsz = d.shape

    # extend r & d to cover entire image. Two possible methods:
    if (always == True):
        # (1) this method lets boundry drift with actual displacements at
        #     edges of 'd' table.

        ns = dsz[1]
        nt = dsz[2]
        Rx, Px = extend(r[0, :, :], d[0, :, :])
        Ry, Py = extend(r[1, :, :], d[1, :, :])

        Ry = np.transpose(Ry)
        Py = np.transpose(Py)

    Ms = np.array([-1,3,-3,1, 3,-6,0,4, -3,3,3,1, 1,0,0,0],
                  order="F")/6.
    Ms = np.reshape(Ms, (4, 4))
    MsT = (Ms)

    sz = scene.shape
    nx = sz[0]
    ny = sz[1]

    ans = np.zeros((nx, ny, 2), order="F")
    for v in range(0, dsz[2]+3):
        t0 = Ry[1, v+1]
        tn = Ry[1, v+2]
        if ((tn < 0) or (t0 > ny-1)):
            break
        t0 = int(np.amax([t0, 0]))
        tn = int(np.amin([tn, ny-1]))
        t = np.arange(tn-t0)/dt + (t0-Ry[1, v+1])/dt
        for u in range(0,dsz[1]+3):
            s0 = Rx[u+1, v+1]
            sn = Rx[u+2, v+1]

            if (sn < 0) or (s0 > nx-1):
                break
            s0 = int(np.amax([s0, 0]))
            sn = int(np.amin([sn, nx-1]))
            s = np.arange(sn-s0)/ds + (s0-Rx[u+1, v+1])/ds
            compx = np.reshape(np.matmul(np.matmul(Ms,
                                                   Px[u:u+4,v:v+4]),
                                         MsT), (4, 4))

            compy = np.reshape(np.matmul(np.matmul(Ms,
                                                   Py[u:u+4,v:v+4]),
                                         MsT), (4, 4))
            ans[s0:sn, t0:tn, :] = patch(compx, compy, s, t)

def patch(compx, compy, s, t):
    """
    TBD

    Parameters
    ----------
    compx : TYPE
        DESCRIPTION.
    compy : TYPE
        DESCRIPTION.
    s : TYPE
        DESCRIPTION.
    t : TYPE
        DESCRIPTION.

    Returns
    -------
    ans: TYPE
        Description

    """
    s     = np.array(s, order="F")
    t     = np.array(t, order="F")

    len_s = len(s)
    len_t = len(t)
    ans = np.zeros((len_s, len_t, 2), order="F")

    ss = np.concatenate((s**3, s**2, s, s**0))
    ss = np.reshape(ss, (len_s, 4))
    tt = np.concatenate((t**3, t**2, t, t**0))
    tt = np.reshape(tt, (len_t, 4)).transpose()

    ans[:, :, 0] = np.matmul(np.matmul(ss, compx), tt)
    ans[:, :, 1] = np.matmul(np.matmul(ss, compy), tt)

    return ans


def extend(cntrlpts_ref, cntrlpts_actl, num_extend_pts=3):
    """Extend map of measured control points and displacements.

    Extend the maps of reference and actual control points to cover area
    outside of measured area. This is necessary to create a smooth displacement
    surface covering the whole scene

    Parameters
    ----------
    cntrlpts_ref : TYPE
        reference control points
    cntrlpts_actl : TYPE
        actual displaced position of control points

    Returns
    -------
    cntrlpts_ref_extnd : TYPE
        reference control points, extended by the appropriate border
    cntrlpts_actl_extnd : TYPE
        actual displaced position of control points, also extended

    """
    # if true, the extended border of the actual displacements will be filled
    #     with the same displacement values as at the corresponding edge.
    # if false, displacements in the extended border will be set to zero,
    #     but that may cause discontinuities in the displacement surface.
    extend_with_offsets = True

    # set the number of additional control points around the edge of the array
    # by which to extend the displacement map
    # num_extend_pts = 3

    # define the size of the extended arrays to generate
    dsz     = cntrlpts_actl.shape
    disp_nx = dsz[0] + num_extend_pts * 2
    disp_ny = dsz[1] + num_extend_pts * 2

    # First, create the entended array of control points
    cntrlpts_ref_extnd = np.zeros((disp_nx, disp_ny), order="F")

    # as currently written, one coordinate of the control point locations are
    # passed into this routine at a time. So the either the values will be varying
    # in the x-direction or the y-direction
    # We compare the differences of the change in values in the two directions
    # to identify which of the two cases we have
    step_x = cntrlpts_ref[1, 0] - cntrlpts_ref[0, 0]
    step_y = cntrlpts_ref[0, 1] - cntrlpts_ref[0, 0]

    # generally, the input arrays will contain the x- or y- coordinates of the control
    # points, so the values will only be varying in one direction if those are laid out
    # on a rectangular grid. In the other direction, the increments between points will 
    # be zero. So we just look to see in which direction is the increment bigger and 
    # use that as the step size to use for the extension.
    # But really it would be better to have the direction to use to define the step size 
    # be defined as an input, or have the step size be an input parameter.
    # This might be useful if the step size is not constant, or changes in both directions
    # simulataneously
    
    # if step_y is greater and non-zero, we'll use that to fill the rows
    if step_x > step_y:
        # define the starting point, which is num_extpts times the step size before the input position
        start_pt    = cntrlpts_ref[0, 0] - num_extend_pts * step_x
        # generate a new array of evenly spaced control points
        new_steps   = np.arange(disp_nx) * step_x + start_pt
        # replicate that array of control points into all rows of the control point array
        for i in range(disp_ny):
            cntrlpts_ref_extnd[:, i] = new_steps
    # if step_y is greater and non-zero, we'll use that to fill the columns
    else:
        # define the starting point, which is num_extpts times the step size before the input position
        start_pt    = cntrlpts_ref[0, 0] - num_extend_pts * step_y
        # generate a new array of evenly spaced control points
        new_steps   = np.arange(disp_ny) * step_y + start_pt
        # replicate that array of control points into all rows of the control point array
        for i in range(disp_nx):
            cntrlpts_ref_extnd[i, :] = new_steps

    # Next create an extended array of the displaced positions of the control points
    # and populate the center portion with the measured (input) offsets

    cntrlpts_actl_extnd = np.zeros((disp_nx, disp_ny), order="F")
    cntrlpts_actl_extnd[num_extend_pts:disp_nx-num_extend_pts, num_extend_pts:disp_ny-num_extend_pts] = cntrlpts_actl - cntrlpts_ref

    # if requested, replicate the edges of the displacement array into the extended boundaries
    if extend_with_offsets is True:
        # extend displacement array by replicating the values of the displacements along the
        #   borders of measured array. This ensures some consistencies in the values

        # take the bottom row of measured offset and replicate it into the extended array of offsets below
        # note: this could be done without a for loop...
        x = cntrlpts_actl_extnd[:, num_extend_pts]
        for i in range(num_extend_pts):
            cntrlpts_actl_extnd[:, i] = x

        # take the top row of measured offset and replicate it into the extended array of offsets above
        x = cntrlpts_actl_extnd[:, disp_ny - num_extend_pts - 1]
        for i in range(disp_ny - num_extend_pts, disp_ny):
            cntrlpts_actl_extnd[:, i] = x

        # take the left column of measured offset and replicate it into the extended array of offsets to the left
        x = cntrlpts_actl_extnd[num_extend_pts, :]
        for i in range(num_extend_pts):
            cntrlpts_actl_extnd[i, :] = x

        # take the left column of measured offset and replicate it into the extended array of offsets to the left
        x = cntrlpts_actl_extnd[disp_nx - num_extend_pts - 1, :]
        for i in range(disp_nx - num_extend_pts, disp_ny):
            cntrlpts_actl_extnd[i, :] = x
        print(cntrlpts_actl_extnd[2, :])

    # now add the control point positions back into the displacement array
    cntrlpts_actl_extnd = cntrlpts_actl_extnd + cntrlpts_ref_extnd

    return cntrlpts_ref_extnd, cntrlpts_actl_extnd


def apod_mask(nx, ny, fraction=0.08):
    """
    Create an apodization mask over the apertures
    to reduce FFT edge effects.

    Parameters
    ----------
    nx : int
        Width of window in pixels
    ny : int
        Height of window in pixels
    fraction: float
        Fraction of window over which intensity drops to zero

    Returns
    -------
    Apodization window (NumPy array)

    """
    taper_wx = int(nx * min(fraction, 0.5))
    taper_wy = int(ny * min(fraction, 0.5))

    filt_x = signal.windows.blackman(2 * taper_wx)
    filt_y = signal.windows.blackman(2 * taper_wy)

    left = filt_x[:taper_wx]
    right = left[::-1]
    top = filt_y[:taper_wy]
    bottom = top[::-1]
    center_x = np.ones(nx - 2*taper_wx)
    center_y = np.ones(ny - 2*taper_wy)

    x_arr = np.concatenate((left, center_x, right))
    y_arr = np.concatenate((top, center_y, bottom))

    m = np.array(np.outer(x_arr, y_arr), order='F')

    return m

def smouth(nx, ny):
    """
    Smouthing window to be applied to the 2D FFTs to
    remove HF noise.

    WORKS! Checked against IDl

    Parameters
    ----------
    nx : integer
        Window size in x-direction.
    ny : integer
        Window size in y-direction.

    Returns
    -------
    mm : ndarry [nx, ny]
        smoothing mask.

    """

    x = np.arange(nx//2)
    if nx % 2 == 1:
        x = np.concatenate([x, x[nx//2-1:nx//2], np.flip(x)])
    else:
        x = np.array([x, np.flip(x)],).flatten()
    if nx > 60:
        magic_number = nx//6
    else:
        magic_number = 10
    x = np.exp(-1*(x/(magic_number))**2)

    y = np.arange(ny//2)
    if (ny % 2) == 1:
        y = np.concatenate([y, y[ny//2-1:ny//2], np.flip(y)])
    else:
        y = np.array([y, np.flip(y)]).flatten()
    if ny > 60:
        magic_number = ny//6
    else:
        magic_number = 10
    y = np.exp(-1*(y/(magic_number))**2)

    mm = np.outer(x.T, y)

    return mm

def doref(ref, apod_mask, d_info, use_fft=False):
    """
    Setup reference window

    Parameters
    ----------
    ref : array (*, *)
        reference image
    mask : TYPE
        mask
    d_info : TYPE
        Destretch_info.

    Returns
    -------
    win: array (*, *)
        Reorganized window
    """

    win = np.zeros((d_info.wx, d_info.wy, d_info.cpx, d_info.cpy),
                   dtype="complex", order="F")
    nelz = d_info.wx * d_info.wy
    ly = d_info.by
    hy = ly + d_info.wy
    for j in range(0, d_info.cpy):
        lx = d_info.bx
        hx = lx + d_info.wx
        for i in range(0, d_info.cpx):
            z = ref[lx:hx, ly:hy]
            z = z - np.sum(z)/nelz
            #z=z-np.polyfit(z[0,  :], z[1:, ],1)
            if use_fft:
                win[:, :, i, j] = np.conj(np.fft.fft2(z*apod_mask))
            else:
                win[:, :, i, j] = z

            lx = lx + d_info.kx
            hx = hx + d_info.kx
        ly = ly + d_info.ky
        hy = hy + d_info.ky


    return win

def cploc(s, w, apod_mask, smou, d_info, use_fft=False, adf2_pad=0.25):
    """
    Locate control points

    Parameters
    ----------
    s : TYPE
        Scene to be registered
    w : TYPE
        Reference image (window) from doref
    mask : TYPE
        DESCRIPTION.
    smou : TYPE
        DESCRIPTION.
    d_info : TYPE
        Destretch information

    Returns
    -------
    ans: TYPE

    """

    ans = np.zeros((2, d_info.cpx, d_info.cpy), order="F")

    nels = d_info.wx * d_info.wy

#    pad_x = int(d_info.kx * adf2_pad)
#    pad_y = int(d_info.ky * adf2_pad)
    pad_x, pad_y = 2, 2

    ly = d_info.by
    hy = ly + d_info.wy
    for j in range(0, d_info.cpy):
        lx = d_info.bx
        hx = lx + d_info.wx

        for i in range(0, d_info.cpx):
            if use_fft:
                #cross correlation, inline
                ss = s[lx:hx, ly:hy]
                #ss = (ss - np.polyfit(ss[0, :], ss[1 1))*mask
                ss_fft = np.array(np.fft.fft2(ss), order="F")
                ss_fft = ss_fft  * w[:, :, i, j] #* smou

                ss_ifft = np.abs(np.fft.ifft2(ss_fft), order="F")
                cc = np.fft.fftshift(ss_ifft)
                #cc = np.roll(ss_ifft, (d_info.wx//2, d_info.wy//2),
                #             axis=(0, 1))

                cc = np.array(cc, order="F")
            else:
                ss = s[lx-pad_x:hx+pad_x, ly-pad_y:hy+pad_y]
                cc = np.zeros((2*pad_x + 1, 2*pad_y + 1))
                for m in range(2*pad_x + 1):
                    for n in range(2*pad_y + 1):
                        cc[m, n] = -np.sum(np.abs(ss[m:m+d_info.wx, n:n+d_info.wy]
                                                 - w[:, :, i, j]))**2
#                cc4 = np.zeros((2*pad_x + 1, 2*pad_y + 1, d_info.wx, d_info.wy))
#                for m in range(2*pad_x + 1):
#                    for n in range(2*pad_y + 1):
#                        cc4[m, n] = ss[m:m+d_info.wx, n:n+d_info.wy]
#                cc = -np.sum(np.abs(cc4 - w[:, :, i, j]), (2, 3))**2
            mx  = np.amax(cc)
            loc = cc.argmax()


            ccsz = cc.shape
            ymax = loc % ccsz[0]
            xmax = loc // ccsz[0]
            # breakpoint()
            #a more complicated interpolation
            #(from Niblack, W: An Introduction to Digital Image Processing, p 139.)

            if ((xmax*ymax > 0) and (xmax < (ccsz[0]-1))
                and (ymax < (ccsz[1]-1))):
            #if (1 == 1):
                denom = 2 * mx - cc[xmax-1,ymax] - cc[xmax+1,ymax]
                xfra = (xmax-1/2) + (mx-cc[xmax-1,ymax])/denom

                denom = 2 * mx - cc[xmax,ymax-1] - cc[xmax,ymax+1]
                yfra = (ymax-1/2) + (mx-cc[xmax,ymax-1])/denom

                # breakpoint()
                xmax=yfra
                ymax=xfra

            if use_fft:
                ans[0,i,j] = lx + xmax
                ans[1,i,j] = ly + ymax
            else:
                ans[0,i,j] = lx + d_info.kx + xmax - pad_x
                ans[1,i,j] = ly + d_info.ky + ymax - pad_y

            lx = lx + d_info.kx
            hx = hx + d_info.kx

        ly = ly + d_info.ky
        hy = hy + d_info.ky


    return ans

def val(scene, ref, kernel):
     # check parameters are reasonable

     # scene, ref, kernel: as defined by 'reg'

     ssz = scene.shape
     rsz = ref.shape
     ksz = kernel.shape
     errflg = 0

     if ((len(ssz) != 2) and (len(ssz) != 3)):
        print("argument 'scene' must be 2-D or 3-D")
        errflg = errflg + 1


     if len(rsz) != 2:
        print("argument 'ref' must be 2-D")
        errflg = errflg + 1

     if ((ssz[0] != rsz[0]) or (ssz[1] != rsz[2])):
         print("arguments 'scene' & 'ref' 1st 2 dimensions must agree")
         errflg = errflg + 1

     if len(ksz) != 2:
        print, "argument kernel must be 2-D"
        errflg = errflg + 1

     if errflg > 0:
         print("quitting - too many errors")

def doreg(scene, r, d, d_info):
    """


    Parameters
    ----------
    scene : TYPE
        Scene to be destretched
    r : TYPE
        reference displacements of the control points
    d : TYPE
        Actual displacements of the control points
    d_info: Destr class
        Destretch information

    Returns
    -------
    ans : TYPE
        Destretched scene.

    """

    xy  = bilin_control_points(scene, r, d)
    ans = bilin_values_scene(scene, xy, d_info)

    return ans

def measure_destr_properties(scene1, scene2):
    """
    Measure the suitable parameters for the destretch based on the two provided
    images based on:
        1)Fourier transforms of the images to obtain the smallest
    scale (box size) on the image;
        2) Fourier transform of the image ratio (???) to measure the
    control point spacing;

    Input:
        -- scene1 -- ndarray (nx, ny)
            Image 1
        -- scene 2 -- ndarray (nx, ny)
            Image 2

    Output:
        -- d_info -- Destretch_class
            Suggested properties of the destretch
    """


    return d_info

def mkcps_nonuniform(ref, kernel):

    return 0

def mkcps_overlapping(ref, kernel, box_size):
    """
    Create control point locations in the reference with overlapping cross
    correlation regions
    """
    return d_info, rcps

def mkcps(ref, kernel, mf=0.08):
    """
    Seems to work
    Choose control point locations in the reference

    Parameters
    ----------
    ref : TYPE
        Reference scene
    kernel : TYPE
        Kernel props

    Returns
    -------
    d_info: Destr class
        Destructor info

    rcps : TYPE
        DESCRIPTION.

    """
    ksz = kernel.shape
    kx = ksz[0]
    ky = ksz[1]


    a = 20./np.sqrt(2.)
    b = 20./np.log(2.)

    wx = int(kx*2)
    wy = int(ky*2)

    if (wx % 2):
        wx = int(wx + 1)
    if (wy % 2):
        wy = int(wy + 1)

    rsz = ref.shape
    cpx = int((rsz[0] - wx + kx)//kx)
    cpy = int((rsz[1] - wy + ky)//ky)

    bx = int(((rsz[0] - wx + kx) % kx)/2)
    by = int(((rsz[1] - wy + ky) % ky)/2)
    rcps = np.zeros((2, cpx, cpy), order="F")
    ly = by
    hy = ly + wy
    for j in range(0, cpy):
        lx = bx
        hx = lx + wx
        for i in range(0, cpx):
            rcps[0, i, j] = (lx + hx)/2
            rcps[1, i, j] = (ly + hy)/2
            lx = lx + kx
            hx = hx + kx

        ly = ly + ky
        hy = hy + ky

    d_info = Destretch_params(kx, ky, wx, wy, bx, by, cpx, cpy, mf)

    return d_info, rcps

def setup(scene, ref, kernel, d_info):

    return

def undo():
    return

def repair(ref, disp, d_info):
    """
    Check if the displacements are good

    Parameters
    ----------
    ref : TYPE
        reference coordinates
    disp : TYPE
        displacements to be checked
    d_info : TYPE
        Destr info

    Returns
    -------
    good : TYPE
        DESCRIPTION.

    """
    TOOFAR = .5             # user may want to change this parameter
    # TOOFAR = 1.0           # user may want to change this parameter
    # TOOFAR = 1.5           # user may want to change this parameter

    sz = disp.shape
    nx = sz[1]
    ny = sz[2]

    if (len(sz) == 4):
        nf = sz[3]
    else:
        nf = 1

    kkx = ref[0, 1, 0] - ref[0, 0, 0]
    kky = ref[1, 0, 1] - ref[1, 0, 0]
    limit = (np.amax([kkx,kky])*TOOFAR)**2

    good = disp + 0

    kps = np.reshape(disp[:, :, :], (2, nx, ny))

    diff = kps - ref

    # list of bad coordinates in this frame
    bad = np.where((diff[0, :, :]**2 + diff[1, :, :]**2) > limit,
                   1, 0)
    bad_count = np.sum(bad)
    i = 0
    j = 0
    while i < bad_count or j < nx*ny:
        x = i % nx
        y = i // nx

        if bad[x, y] == 1:
            good [:, x, y] = ref[:, x, y]
            i += 1
        j += 1
    return good

def cps(scene, ref, kernel, use_fft=False, adf2_pad=0.25):
    """
    Control points for sequence destretch

    Parameters
    ----------
    scene : TYPE
        Input scene [nx, ny] or [nx, ny, nf] for which
        displacements are computed
    ref : TYPE
        Reference scene [nx, ny]
    kernel : TYPE
        Kernel size [kx, ky]
    d_info : TYPE
        DESCRIPTION.

    Returns
    -------
    ans : TYPE
        displacement array [2, cpx, cpy, nf]

    """

    d_info, rdisp = mkcps(ref, kernel)


    #mm = np.zeros((d_info.wx,d_info.wy), order="F")
    #mm[:, :] = 1
    mm = apod_mask(d_info.wx, d_info.wy, d_info.mf)

    smou = np.zeros((d_info.wx,d_info.wy), order="F")
    smou[:, :] = 1


    ref = ref/np.average(ref)*np.average(scene)
    win = doref(ref, mm, d_info, use_fft)

    ssz = scene.shape
    nf = ssz[2]
    ans = np.zeros((2, d_info.cpx, d_info.cpy, nf), order="F")



    # compute control point locations
    for frm in range(0, nf):
        ans[:, :, :, frm] = cploc(scene[:, :, :, frm], win, mm, smou, d_info, use_fft, adf2_pad)
        #ans[:, :, :, frm] = repair(rdisp, ans[:, :, :,frm], d_info)# optional repair

    if ssz[2]:
        scene = np.reshape(scene, (ssz[0], ssz[1]))
        ans = np.reshape(ans, (2, d_info.cpx, d_info.cpy))


    return ans

def reg(scene, ref, kernel_size, mf=0.08, use_fft=False, adf2_pad=0.25):
    """
    Register scenes with respect to ref using kernel size and
    then returns the destretched scene.

    Parameters
    ----------
    scene : [nx, ny] [nx, ny, nf]
        Scene to be registered
    ref : [nx, ny]
        reference frame
    kernel_size : int
       Kernel size (otherwise unused)!!!!!

    Returns
    -------
    ans : [nx, ny]
        Destreched scene.
    disp : ndarray (kx, ky)
        Control point locations
    rdisp : ndarray (kx, ky)
        Reference control point locations

    """
    scene -= scene.mean()
    ref -= ref.mean()
    kernel = np.zeros((kernel_size, kernel_size))

    d_info, rdisp = mkcps(ref, kernel, mf)
    mm = apod_mask(d_info.wx, d_info.wy, d_info.mf)
    smou = smouth(d_info.wx, d_info.wy)
    #Condition the ref
    win = doref(ref, mm, d_info, use_fft)


    ssz = scene.shape
    ans = np.zeros((ssz[0], ssz[1]), order="F")

    # compute control point locations

    #start = time()
    disp = cploc(scene, win, mm, smou, d_info, use_fft, adf2_pad)
    #end = time()
    #dtime = end - start
    #print(f"Time for a scene destretch is {dtime:.3f}")

    #disp = repair(rdisp, disp, d_info) # optional repair
    #rms = sqrt(total((rdisp - disp)^2)/n_elements(rdisp))
    #print, 'rms =', rms
    #mdisp = np.mean(rdisp-disp,axis=(1, 2))
    #disp[0, :, :] += mdisp[0]
    #disp[1, :, :] += mdisp[1]
    x = doreg(scene, rdisp, disp, d_info)
    ans = x
        #    win = doref (x, mm); optional update of window

#    print(f"Total destr took: {(end - start):.5f} seconds for kernel"
 #         +f"of size {kernel_size} px.")

    return ans, disp, rdisp, d_info

def reg_saved_window(scene, win, kernel_size, d_info, rdisp, mm, smou, use_fft=False, adf2_pad=0.25):
    """
    Register scenes with respect to ref using kernel size and
    then returns the destretched scene, using precomputed window.

    Parameters
    ----------
    scene : [nx, ny] [nx, ny, nf]
        Scene to be registered
    win: [nx, ny, nf]
        FFT of the reference scene (computed with doref)
    kernel_size : int
       Kernel size (otherwise unused)!!!!!

    Returns
    -------
    ans : [nx, ny]
        Destreched scene.
    disp : ndarray (kx, ky)
        Control point locations
    rdisp : ndarray (kx, ky)
        Reference control point locations

    """
    kernel = np.zeros((kernel_size, kernel_size))

    # d_info, rdisp = mkcps(ref, kernel)
    # mm = mask(d_info.wx, d_info.wy)
    # smou = smouth(d_info.wx, d_info.wy)
    #Condition the ref
    #win = doref(ref, mm, d_info)


    ssz = scene.shape
    ans = np.zeros((ssz[0], ssz[1]), order="F")

    # compute control point locations

    #start = time()
    disp = cploc(scene, win, mm, smou, d_info, use_fft, adf2_pad)
   # end = time()
  #  dtime = end - start
 #   print(f"Time for a scene destretch is {dtime:.3f}")

    #disp = repair(rdisp, disp, d_info) # optional repair
    #rms = sqrt(total((rdisp - disp)^2)/n_elements(rdisp))
    #print, 'rms =', rms
    #mdisp = np.mean(rdisp-disp,axis=(1, 2))
    #disp[0, :, :] += mdisp[0]
    #disp[1, :, :] += mdisp[1]
    x = doreg(scene, rdisp, disp, d_info)
    ans = x
        #    win = doref (x, mm); optional update of window




    #print(f"Total destr took: {(end - start):.5f} seconds for kernel"
    #      +f"of size {kernel_size} px.")

    return ans, disp, rdisp, d_info


def reg_loop(scene, ref, kernel_sizes, mf=0.08, use_fft=False, adf2_pad=0.25):
    """
    Parameters
    ----------
    scene : ndarray (nx, ny)
        Image to be destretched
    ref : ndarray (nx, ny)
        Reference image
    kernel_sizes : ndarray (n_kernels)
        Sizes of the consecutive kernels to be applied

    Returns
    -------
    ans : ndarray (nx, ny)
        Destretched scene
    d_info: Destretch class
        Parameters of the destretching
    """


    scene_temp = scene
    start = time()

    for el in kernel_sizes:
        scene_temp, disp, rdisp, d_info = reg(scene_temp, ref, el, mf, use_fft, adf2_pad)

    end = time()
    print(f"Total elapsed time {(end - start):.4f} seconds.")
    ans = scene_temp

    return ans, disp, rdisp, d_info


def reg_loop_series(scene, ref, kernel_sizes, mf=0.08, use_fft=False, adf2_pad=0.25):
    """
    Parameters
    ----------
    scene : ndarray (nx, ny, nt)
        Image to be destretched
    ref : ndarray (nx, ny)
        Reference image
    kernel_sizes : ndarray (n_kernels)
        Sizes of the consecutive kernels to be applied

    Returns
    -------
    ans : ndarray (nx, ny)
        Destretched scene
    d_info: Destretch class
        Parameters of the destretching
    """

    num_scenes = scene.shape[2]
    scene_d = np.zeros((scene.shape))

    start = time()
    num_kernels = len(kernel_sizes)
    windows = {}
    d_info_d = {}
    mm_d = {}
    smou_d = {}
    rdisp_d = {}

    # d_info, rdisp = mkcps(ref, kernel)
    # mm = mask(d_info.wx, d_info.wy)
    # smou = smouth(d_info.wx, d_info.wy)
    for kernel1 in kernel_sizes:
        kernel = np.zeros((kernel1, kernel1))

        d_info, rdisp = mkcps(ref, kernel, mf)
        d_info_d[kernel1] = d_info
        rdisp_d[kernel1] = rdisp

        mm = apod_mask(d_info.wx, d_info.wy, d_info.mf)
        mm_d[kernel1] = mm

        smou = smouth(d_info.wx, d_info.wy)
        smou_d[kernel1] = smou

        win = doref(ref, mm, d_info, use_fft)
        windows[kernel1] = win
        
    disp_l = list(rdisp.shape)
    disp_l.append(num_scenes)
    disp_t = tuple(disp_l)
    disp_all = np.zeros(disp_t)

    for t in range(num_scenes):
        for k in kernel_sizes:
            scene_d[:, :, t], disp, rdisp, d_info = reg_saved_window(scene[:, :, t],
                                                                     windows[k],
                                                                     k, d_info_d[k],
                                                                     rdisp_d[k],
                                                                     mm_d[k],
                                                                     smou_d[k],
                                                                     use_fft,
                                                                     adf2_pad)
        disp_all[:, :, :, t] = disp

    end = time()
    print(f"Total elapsed time {(end - start):.4f} seconds.")
    ans = scene_d

    return ans, disp_all, rdisp, d_info

def test_destretch(scene, ref, kernel_size, plot=False):
    start = time()
    ans1, disp, rdisp, d_info = reg_loop(scene, ref, kernel_size)
    if plot==True:
        pl.figure(dpi=250)
        pl.imshow(scene, origin=0)
        pl.title("Original scene")
        pl.show()

        pl.figure(dpi=250)
        pl.imshow(ans1, origin=0)
        pl.title("Destretched scene")
        pl.show()

        pl.figure(dpi=250)
        pl.imshow(ref, origin=0)
        pl.title("Reference")
        pl.show()
    end = time()
    print(f"Total elapsed time for test_function is {end-start}.")

def test_rotation(scene, angle):
    """
        Test if the code can pick up a static rotation of an image
    """



def test_residual_diff():
    """
        Test if the code reduces the distortions between the two images
    """
