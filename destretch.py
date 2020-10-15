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

class Destretch_params():
    
    def __init__(self, kx, ky, wx, wy, bx, by, cpx, cpy):
        self.kx = kx
        self.ky = ky
        self.wx = wx
        self.wy = wy
        self.bx = bx
        self.by = by
        self.cpx = cpx
        self.cpy = cpy
    
    def print_props(self):
        print("[kx, ky, wx, wy, bx, by, cpx, cpy] are:")
        print(self.kx, self.ky, self.wx, self.wy, self.bx, self.by,
              self.cpx, self.cpy)
              
def bilin(s, xy, d_info, nearest_neighbor = False):
    """
    Bilinear interpolation of the scene s

    Parameters
    ----------
    s : TYPE
        Scene
    xy : TYPE
        ???
    d_info: class Destretch_Params
        Destretch parameters

    Returns
    -------
    ans: 
        Bilinear interpolated 

    """
    
    if nearest_neighbor == True:
        x = np.array(xy[:, :, 0] + .5, order="F")
        y = np.array(xy[:, :, 1] + .5, order="F")    
    
    else:
        x = np.array(xy[:, :, 0], order="F")
        y = np.array(xy[:, :, 1], order="F")
        
        x0 = x.astype(int)
        x1 = (x+1).astype(int)
        y0 = (y).astype(int)
        y1 = (y+1).astype(int)
        
        fx = x % 1.
        fy = y % 1.
        
        s  = np.array(s, order="F")
        ss = s.astype(float)
        
        ss00 = ss[x0, y0]
        ss01 = ss[x0, y1]
        ssfx = (ss[x1, y0] - ss00) * fx
        ssfy = fy * (ss01 - ss00 + (ss[x1, y1] - ss01) * fx - ssfx)
        ans  = ss00 + ssfx + ssfy
    
    return ans

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
        print(t)
        for u in range(0,dsz[1]+3):
            s0 = Rx[u+1, v+1]
            sn = Rx[u+2, v+1]
            
            if (sn < 0) or (s0 > nx-1):
                break 
            s0 = int(np.amax([s0, 0]))
            sn = int(np.amin([sn, nx-1]))
            s = np.arange(sn-s0)/ds + (s0-Rx[u+1, v+1])/ds
            print(s0, sn, t0, tn)
            compx = np.reshape(np.matmul(np.matmul(Ms, 
                                                   Px[u:u+4,v:v+4]),
                                         MsT), (4, 4))
            print(compx.shape)
            compy = np.reshape(np.matmul(np.matmul(Ms, 
                                                   Py[u:u+4,v:v+4]),
                                         MsT), (4, 4))
            ans[s0:sn, t0:tn, :] = patch(compx, compy, s, t)

    return ans

def mask(nx, ny):
    """
    Create a mask over the apertures

    Parameters
    ----------
    nx : TYPE
        DESCRIPTION.
    ny : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    m = np.ones((int(nx), int(ny)), order="F")
    
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
    hy = ly + d_info.wy -1 
    for j in range(0, d_info.cpy):
        lx = d_info.bx
        hx = lx + d_info.wx - 1
        for i in range(0, d_info.cpx):
            z = ref[lx:hx+1, ly:hy+1]
            z = z - np.sum(z)/nelz
            #z=z-np.polyfit(z[0,  :], z[1:, ],1)
            win[:, :, i, j] = np.conj(np.fft.fft2(z*mask))

            lx = lx + d_info.kx
            hx = hx + d_info.kx

        ly = ly + d_info.ky
        hy = hy + d_info.ky

    return win
    
def cploc(s, w, mask, smou, d_info):
    """
    Locate control points

    Parameters
    ----------
    s : TYPE
        Scene to be registered
    w : TYPE
        Reference image from doref
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
    
    ly = d_info.by
    hy = ly + d_info.wy
    for j in range(0, d_info.cpy):
        lx = d_info.bx
        hx = lx + d_info.wx
        
        for i in range(0, d_info.cpx):

            #cross correlation, inline
            ss = s[lx:hx, ly:hy]

            #ss = (ss - np.polyfit(ss[0, :], ss[1 1))*mask
            ss_fft = np.fft.fft2(ss) * w[:, :, i, j] * smou
            ss_ifft = np.abs(np.fft.ifft2(ss_fft))
            cc = np.roll(ss_ifft, (d_info.wx//2, d_info.wy//2), 
                         axis=(0, 1))

            cc = np.array(cc, order="F")
            mx  = np.amax(cc)
            loc = cc.argmax()
            
        
            ccsz = cc.shape
            xmax = loc % ccsz[0]
            ymax = loc // ccsz[0]

            #a more complicated interpolation
            #(from Niblack, W: An Introduction to Digital Image Processing, p 139.)
            
            if ((xmax*ymax > 0) and (xmax < (ccsz[0]-1)) 
                and (ymax < (ccsz[1]-1))):
    
                denom = mx*2 - cc[xmax-1,ymax] - cc[xmax+1,ymax]
                xfra = (xmax-.5) + (mx-cc[xmax-1,ymax])/denom

                denom = mx*2 - cc[xmax,ymax-1] - cc[xmax,ymax+1]
                yfra = (ymax-.5) + (mx-cc[xmax,ymax-1])/denom

                xmax=xfra
                ymax=yfra
            

            ans[0,i,j] = lx + xmax
            ans[1,i,j] = ly + ymax
        
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
    
    xy = bspline(scene, r, d, d_info,)
    ans = bilin(scene, xy, d_info)
    
    return ans
         
def mkcps(ref, kernel):
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
    
    d_info = Destretch_params(kx, ky, wx, wy, bx, by, cpx, cpy) 
    
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
    print(disp.shape)

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


    mm=np.zeros((d_info.wx,d_info.wy), order="F")
    mm[:, :] = 1

    smou=np.zeros((d_info.wx,d_info.wy), order="F")
    smou[:, :] = 1


    ref = ref/np.average(ref)*np.average(scene)
    win = doref(ref, mm, d_info)

    ssz = scene.shape
    nf = ssz[2]
    ans = np.zeros((2, d_info.cpx, d_info.cpy, nf), order="F")



    # compute control point locations
    for frm in range(0, nf):
        ans[:, :, :, frm] = cploc(scene[:, :, :, frm], win, mm, smou, d_info)
        ans[:, :, :, frm] = repair(rdisp, ans[:, :, :,frm], d_info)# optional repair
    
    if ssz[2]:
        scene = np.reshape(scene, (ssz[0], ssz[1]))
        ans = np.reshape(ans, (2, d_info.cpx, d_info.cpy))
    
    return ans

def reg(scene, ref, kernel, rdisp, disp, apply_scene=0):
    """
    Register scenes with respect to ref using kernel size and
    then returns the destretched scene.

    Parameters
    ----------
    scene : [nx, ny] [nx, ny, nf]
        Scene to be registered
    ref : [nx, ny]
        reference frame
    kernel : [kx, ky] bit array == np.zeros((kx, ky))
       Kernel size (otherwise unused)!!!!!
    apply_scene : TYPE
        DESCRIPTION.
    rdisp : TYPE
        DESCRIPTION.
    disp : TYPE
        DESCRIPTION.

    Returns
    -------
    ans : [nx, ny, nf]
        Destreched scene.

    """
    
    if apply_scene == 0:
        apply_scene = scene
        
    d_info, rdisp = mkcps(ref, kernel)
    mm = mask(d_info.wx, d_info.wy)
    smou = smouth(d_info.wx, d_info.wy)

    #Condition the ref
    win = doref(ref, mm, d_info)

    ssz = scene.shape
    ans = np.zeros((ssz[0], ssz[1]), order="F")

    # compute control point locations

    disp = cploc(scene[:, :], win, mm, smou, d_info)
    #disp = repair(rdisp, disp, d_info) # optional repair
    #rms = sqrt(total((rdisp - disp)^2)/n_elements(rdisp))
    #print, 'rms =', rms
    x = doreg(apply_scene[:, :], rdisp, disp, d_info)
    ans[:, :] = x
        #    win = doref (x, mm); optional update of window

    undo()

    
    return ans, disp, rdisp



