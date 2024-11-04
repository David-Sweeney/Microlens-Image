import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.constants as cst
import astropy.units as u
from scipy.interpolate import RegularGridInterpolator
from skimage.measure import block_reduce

axes_size = 14
interpolate_value = 0.
    
M_lens = 8*u.M_sun
D_os = 8*u.kpc
D_ol = 4*u.kpc
D_ls = D_os - D_ol

theta_E = np.sqrt((4*cst.G*M_lens)/(cst.c**2) * (D_ls/(D_ol*D_os)))
theta_E = (theta_E * u.rad).to(u.mas)

def load_image(filepath, downsample=None):
    im = plt.imread(filepath)
    if downsample is not None and downsample > 1:
        im = block_reduce(im, (downsample, downsample, 1), np.mean)
    
    if im.max() > 1.:
        im = im/255.
    
    im = np.transpose(im, (1, 0, 2))
    im = np.flip(im, axis=1)
    
    return im

def lens_to_source_plane(theta, mass, D_L, D_LS):
    sep = np.linalg.norm(theta, axis=0).to(u.rad).value
    return theta*(1 - (4*cst.G*mass)/(cst.c**2*D_L * sep**2) * (D_LS/D_L))

def get_axes(im):
    
    xs = np.arange(im.shape[0])
    xs = xs - xs.mean()
    old_x_max = xs.max()
    xs = xs / len(xs) * axes_size
    
    ys = np.arange(im.shape[1])
    ys = ys - ys.mean()
    ys = ys / old_x_max * xs.max()
    
    return xs, ys

def get_interpolators(image, xs, ys, interpolate_value):

    interpolators = []
    for i in range(image.shape[-1]):
        interpolators.append(RegularGridInterpolator((xs, ys), image[..., i], fill_value=interpolate_value, bounds_error=False))

    return interpolators

def interpolate(interpolators, xx, yy):

    output = np.zeros((*xx.shape, 3))

    for i in range(3):
        output[..., i] = interpolators[i]((xx, yy))

    output = np.clip(output, 0., 1.)
    
    return output

def apply_lensing(BH_centre, interpolators, xv, yv):
    # Convert grid to coordinates centred on BH
    xx = xv - BH_centre[0]
    yy = yv - BH_centre[1]
    
    # Calculate lens distortion
    xx, yy = lens_to_source_plane(np.array([xx, yy])*u.mas, M_lens, D_ol, D_ls)
    
    # Convert grid back to image coordinates
    xx = (xx + BH_centre[0]).to(u.mas).value
    yy = (yy + BH_centre[1]).to(u.mas).value
    
    warped = interpolate(interpolators, xx, yy)
    return warped

def warp_image(im, BH_centre, filepath):
    # Train interpolators on original image
    xs, ys = get_axes(im)
    interpolators = get_interpolators(im, xs, ys, interpolate_value=0.)
    extent = np.array([xs.min(), xs.max(), ys.min(), ys.max()])
    
    # Set up figure with no border and correct dimensions
    fig, ax = plt.subplots(figsize=(5, 5*im.shape[1]/im.shape[0]), dpi=512)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.axis('off')
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    
    # Get grid over image
    xv, yv = np.meshgrid(xs, ys)*u.mas

    warped = apply_lensing(BH_centre, interpolators, xv, yv)

    ax.imshow(warped, extent=extent, origin='lower')
    
    multipliers = np.geomspace(1, 1.2, 20)
    alphas = np.exp(-3*(np.arange(len(multipliers))/len(multipliers)))
    for multiplier, alpha in zip(multipliers, alphas):
        ax.add_patch(mpl.patches.Circle(BH_centre.to(u.mas).value, radius=1.9*multiplier, alpha=alpha, fill=True, fc='k', ec=None))
    plt.savefig(filepath)
    
if __name__ == '__main__':
    assert len(sys.argv) == 3, "Usage: python microlens_image.py <image> <output>"
    
    # Load image from command line argument
    im = load_image(sys.argv[1])
    
    warp_image(im, np.array([0, 0])*u.mas, sys.argv[2])
    
