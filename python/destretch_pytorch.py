#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 16:09:55 2020

Destretching routines for removing optical defects following the
reg.pro by Phil Wiborg and Thomas Rimmele in reg.pro

@author: molnarad
"""

import numpy as np
import scipy as sp
from scipy import signal as signal
from scipy.ndimage.interpolation import shift
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as pl
from time import time
import torch as t
import torch.fft as fft
from flicker import flicker 

class Destretch_params():
    """
    Class containing all the information about the 
    destretching procedure.

    """

    def __init__(self, kx, ky, wx, wy, bx, by, cpx, cpy, device):
        "
        Initialize the class with the following parameters:
            - kx - integer --
        "

        """
        Initiate the Destretch class with the following parameters
        :param kx: spacing of the control points in x
        :param ky: spacing of the control points in y
        :param wx:
        :param wy:
        :param bx:
        :param by:
        :param cpx:
        :param cpy:
        :param device: device to be used (CPU or GPU)
        """
        self.kx = kx
        self.ky = ky
        self.wx = wx
        self.wy = wy
        self.bx = bx
        self.by = by
        self.cpx = cpx
        self.cpy = cpy
        self.device = device

    def print_props(self):
        """
        Print the destretching properties
        :return:
        """
        print("[kx, ky, wx, wy, bx, by, cpx, cpy] are:")
        print(self.kx, self.ky, self.wx, self.wy, self.bx, self.by,
              self.cpx, self.cpy)

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

def bilin_values_scene(s, xy, d_info, nearest_neighbor = False):
    """
    Bilinear interpolation (resampling)
    of the scene s at coordinates xy

    Parameters
    ----------
    s : ndarray (nx, ny)
        Scene
    xy : ndarray (2, nx, ny)
        coordinates of the pixels of the output image
        on the input image (at which to interpolate the scene)
    d_info: class Destretch_Params
        Destretch parameters

    Returns
    -------
    ans: ndarray (nx, ny)
        Bilinear interpolated (resampled) image at the xy locations

    """

    if nearest_neighbor == True:
        x = np.array(xy[:, :, 0] + .5, order="F")
        y = np.array(xy[:, :, 1] + .5, order="F")

    else:
        
        dtype_long = t.cuda.LongTensor
        dtype = t.cuda.FloatTensor

        x0 = t.floor(xy[0, :, :]).type(dtype_long)
        x1 = x0 + 1        
        
        y0 = t.floor(xy[1, :, :]).type(dtype_long)
        y1 = y0 + 1 
       
        x0 = t.clamp(x0, 0, s.shape[1]-1)
        x1 = t.clamp(x1, 0, s.shape[1]-1)
        y0 = t.clamp(y0, 0, s.shape[0]-1)
        y1 = t.clamp(y1, 0, s.shape[0]-1)
     
        fx = xy[0, :, :] % 1
        fy = xy[0, :, :] % 1 

        # print(s, x0.shape)
        ss00 = s[x0, y0]
        ss01 = s[x0, y1]
        ssfx = (s[x1, y0] - ss00) * fx
        ssfy = fy * (ss01 - ss00 + (s[x1, y1] - ss01) * fx - ssfx)
        ans  = ss00 + ssfx + ssfy

        
    return ans

def bilin_control_points(scene, rdisp, disp, d_info,
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
 
    xy_ref_coordinates = np.zeros((2, scene_nx, scene_ny), order="F")



    
    xy_ref_coordinates[0, :, :] = [np.linspace(0, (scene_nx-1), 
                                               num=scene_ny, dtype="int") 
                                   for el in range(scene_nx)]
    xy_ref_coordinates[1, :, :] = [np.zeros(scene_ny, dtype="int")+el 
                                   for el in range(scene_nx)]
    
    xy_ref_coordinates = np.swapaxes(xy_ref_coordinates, 1, 2)
   
  
    xy_ref_coordinates = t.tensor(xy_ref_coordinates, 
                                  device=d_info.device)

    dd = disp - rdisp

    cp_x_coords = cp_x_coords.clone().detach().to("cpu")
    cp_y_coords = cp_y_coords.clone().detach().to("cpu")
    
    dd = dd.clone().detach().to("cpu")
 
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

    xy_grid = t.tensor(xy_grid, device=d_info.device) 
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
    s = np.array(s, order="F")
    t = np.array(t, order="F")

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

def extend(r, d):
    """
    Extend reference and actual displacements to cover whole
    scene

    Parameters
    ----------
    r : TYPE
        reference displacement of points
    d : TYPE
        actual displacements of points


    Returns
    -------
    rd : TYPE
        DESCRIPTION.
    sd : TYPE
        DESCRIPTION.

    """

    dsz = d.shape
    ns  = dsz[0] + 6
    nt  = dsz[1] + 6

    rd = np.zeros((ns, nt), order="F")
    dif = r[1, 0] - r[0, 0]
    zro = r[0, 0] - 3*dif
    z   = np.arange(ns)*dif + zro
    for i in range(nt):
        rd[:, i] = z

    sd = np.zeros((ns, nt), order="F")
    sd[3:ns-3, 3:nt-3] = d -r

    x = sd[:, 3]
    sd[:, 0] = x
    sd[:, 1] = x
    sd[:, 2] = x

    x = sd[: ,nt-4]
    sd[:, nt-3] = x
    sd[:, nt-2] = x
    sd[:, nt-1] = x

    sd = np.transpose(sd)

    x = sd[:, 3]
    sd[:, 0] = x
    sd[:, 1] = x
    sd[:, 2] = x

    x = sd[: ,ns-4]
    sd[:, ns-3] = x
    sd[:, ns-2] = x
    sd[:, ns-1] = x

    sd = np.transpose(sd)
    sd = sd + rd

    return rd, sd

def mask(nx, ny, device):
    """
    Create a mask over the apertures

    Parameters
    ----------
    nx : TYPE
        DESCRIPTION.
    ny : TYPE
        DESCRIPTION.
    device: PyTorch Object
        Tells you if you should use the CPU or GPU
    Returns
    -------
    None.

    """
    m = t.ones((int(nx), int(ny)), device=device)

    return m

def smouth(nx, ny):
    """
    WORKS! Checked against IDl

    Parameters
    ----------
    nx : TYPE
        DESCRIPTION.
    ny : TYPE
        DESCRIPTION.

    Returns
    -------
    mm : TYPE
        DESCRIPTION.

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

def doref(ref, mask, d_info):
    """
    Setup reference window
    And load it to the gpu
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
   
    z = t.zeros((d_info.cpx, d_info.cpy, d_info.wx, d_info.wy,), 
                device=d_info.device)
    win1 = t.zeros(d_info.wx, d_info.wy, d_info.cpx, d_info.cpy)
    win = t.complex(win1, win1)
    win = win.to(d_info.device)
    nelz = d_info.wx * d_info.wy
    ly = d_info.by
    hy = ly + d_info.wy -1
    for j in range(0, d_info.cpy):
        lx = d_info.bx
        hx = lx + d_info.wx - 1
        for i in range(0, d_info.cpx):
            z[i, j, :, :] = ref[lx:hx+1, ly:hy+1]
            #z=z-np.polyfit(z[0,  :], z[1:, ],1)
            lx = lx + d_info.kx
            hx = hx + d_info.kx
        ly = ly + d_info.ky
        hy = hy + d_info.ky
    
    z_sum = t.sum(z, dim=(2, 3))/nelz
    z = z - z_sum[:, :, None, None] 
    z = z * mask
    win = t.conj(t.fft.fftn(z, dim=(2, 3)))

    # Win seems to have jumbled dimensions compared to destretch.py
    # as win have dimensions of [cpx, cpy, wbx, wby] 

    return win

def cploc(s, w, mask, smou, d_info):
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
         
    nx_s = s.shape[0]
    ny_s = s.shape[0]
    ans = t.zeros((2, d_info.cpx, d_info.cpy), device=d_info.device)
    
   # s_reshape = t.reshape(s, (d_info.cpx, d_info.cpy, 
   #                          d_info.kx, d_info.ky))     

    s_vector = t.zeros((d_info.cpx, d_info.cpy, 
                        d_info.wx, d_info.wy), device=d_info.device)

    nelz = d_info.wx * d_info.wy
    ly = d_info.by
    hy = ly + d_info.wy -1
    for j in range(0, d_info.cpy):
        lx = d_info.bx
        hx = lx + d_info.wx - 1
        for i in range(0, d_info.cpx):
            s_vector[i, j, :, :] = s[lx:hx+1, ly:hy+1]
            #z=z-np.polyfit(z[0,  :], z[1:, ],1)
            lx = lx + d_info.kx
            hx = hx + d_info.kx
        ly = ly + d_info.ky
        hy = hy + d_info.ky
    

    s_fft = fft.fftn(s_vector, dim=(2, 3))
    ss_fft = s_fft * w * smou
    
    ss_ifft = t.abs(t.fft.ifftn(ss_fft, dim=(2, 3)))

    # print(ss_ifft.shape)
    cc = t.roll(ss_ifft, (d_info.wx//2, d_info.wy//2), dims=(2, 3)) 

            
    mx  = t.amax(cc, dim=(2, 3))
    ss_ifft_flat = t.flatten(cc, start_dim=2, end_dim=3)
    loc          = t.argmax(ss_ifft_flat, dim=2) 
    ccsz = cc.shape
    xmax1 = loc % ccsz[-1]
    ymax1 = loc // ccsz[-2]

    xmax = xmax1.to("cpu") 
    ymax = ymax1.to("cpu")
    
    x_mesh = ((np.mgrid[0:d_info.cpx:1, 0:d_info.cpy:1])[0]*d_info.kx 
              + d_info.bx)
    y_mesh = ((np.mgrid[0:d_info.cpx:1, 0:d_info.cpy:1])[1]*d_info.ky 
              + d_info.bx)
    
    ans[0,:, :] = t.from_numpy(x_mesh) + xmax
    ans[1,:, :] = t.from_numpy(y_mesh) + ymax

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

    xy  = bilin_control_points(scene, r, d, d_info)
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

def mkcps(ref, kernel, device="cpu"):
    """
    Seems to work
    Choose control point locations in the reference

    Parameters
    ----------
    ref : TYPE
        Reference scene
    kernel : TYPE
        Kernel props
    device: PyTorch device object
        Tells you if you're using the cpu or the gpu
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
    rcps = t.zeros((2, cpx, cpy), device=device)
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

    d_info = Destretch_params(kx, ky, wx, wy, bx, by, cpx, cpy,
                              device)

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

def cps(scene, ref, kernel):
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


    mm = np.zeros((d_info.wx,d_info.wy), order="F")
    mm[:, :] = 1

    smou = np.zeros((d_info.wx,d_info.wy), order="F")
    smou[:, :] = 1


    ref = ref/np.average(ref)*np.average(scene)
    win = doref(ref, mm, d_info)

    ssz = scene.shape
    nf = ssz[2]
    ans = np.zeros((2, d_info.cpx, d_info.cpy, nf), order="F")



    # compute control point locations
    for frm in range(0, nf):
        ans[:, :, :, frm] = cploc(scene[:, :, :, frm], win, mm, smou, d_info)
        #ans[:, :, :, frm] = repair(rdisp, ans[:, :, :,frm], d_info)# optional repair

    if ssz[2]:
        scene = np.reshape(scene, (ssz[0], ssz[1]))
        ans = np.reshape(ans, (2, d_info.cpx, d_info.cpy))


    return ans

def reg(scene, ref, kernel_size, device=False):
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

    kernel = t.zeros((kernel_size, kernel_size), device=device)
    
    
    d_info, rdisp = mkcps(ref, kernel, device)

    mm = mask(d_info.wx, d_info.wy, device)
    mm = mm.clone().detach().to(device)
  
    smou = smouth(d_info.wx, d_info.wy)
    smou = t.tensor(smou, device=device) 
   
    #Condition the ref
    win = doref(ref, mm, d_info) 
    win = win.clone().detach().to(device)
    
    ssz = scene.shape
    
    ans = np.zeros((ssz[0], ssz[1]), order="F")
    
    # compute control point locations

    start = time()
    disp = cploc(scene, win, mm, smou, d_info)
    end = time()
    #disp = repair(rdisp, disp, d_info) # optional repair
    #rms = sqrt(total((rdisp - disp)^2)/n_elements(rdisp))
    #print, 'rms =', rms
    #mdisp = np.mean(rdisp-disp,axis=(1, 2))
    #disp[0, :, :] += mdisp[0]
    #disp[1, :, :] += mdisp[1]

    x = doreg(scene, rdisp, disp, d_info)

    ans = x
        #    win = doref (x, mm); optional update of window 

    
    breakpoint()
    print(f"Total destr took: {end - start} seconds for kernel"
          +f"of size {kernel_size} px.")

    return ans, disp, rdisp, d_info

def reg_loop(scene, ref, kernel_sizes, device="False"):
    """


    Parameters
    ----------
    scene : ndarray (nx, ny)
        Image to be destretched
    ref : ndarray (nx, ny)
        Reference image
    kernel_sizes : ndarray (n_kernels)
        Sizes of the consecutive kernels to be applied
    device: string
        == "CPU" for using the cpu 
        else will try to use the GPU
    Returns
    -------
    ans : ndarray (nx, ny)
        Destretched scene
    d_info: Destretch class
        Parameters of the destretching
    """
   
   
    if t.cuda.is_available():  
        dev = "cuda:0" 
        t.backends.cudnn.benchmark = True
    else:  
        dev = "cpu" 
 
    
    if device == "CPU":
        device = t.device("cpu") 
    else: 
        device = t.device(dev) 

    scene = t.tensor(scene, device=device)
    ref   = t.tensor(ref, device=device)
   
    for el in kernel_sizes:
        scene1, disp, rdisp, d_info = reg(scene, ref, el, 
                                         device=device)

   
    ans = scene1.to("cpu")
    scene = scene.to("cpu")
    ref = ref.to("cpu") 
    
    

    return ans, disp, rdisp, d_info

def test_destretch(scene, ref, kernel_size, plot=False, device=False):

    start = time() 
    scene1 = scene.copy()
    ref1   = ref.copy()
    ans1, disp, rdisp, d_info = reg_loop(scene, ref, kernel_size, 
                                         device=device)
    if plot==True:
        pl.figure(dpi=250)
        pl.imshow(scene1, origin=0)
        pl.title("Original scene")
        pl.show()

        pl.figure(dpi=250)
        pl.imshow(ans1, origin=0)
        pl.title("Destretched scene")
        pl.show()

        pl.figure(dpi=250)
        pl.imshow(ref1, origin=0)
        pl.title("Reference")
        pl.show()

        flicker(scene1, ans1)
        flicker(ref1, ans1)
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
