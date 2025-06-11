# Finish writing moveFrame function
import os
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.constants as cst
import astropy.units as u
from scipy.interpolate import RegularGridInterpolator
from skimage.measure import block_reduce
from lensing_event_rate import getProperMotion

# Make warping image an animation: try mpl.animation (matplotlib.animation)
# Move BH_centre across image
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter

axes_size = 14
interpolate_value = 0.

M_lens = 8*u.M_sun # Mass of lensing object  # Greater mass = greater lensing
D_os = 1000*u.kpc # Distance to background source # Greater distance = greater lensing
D_ol = 80*u.pc # Distance to lensing object # Greater distnace = less lensing
D_ls = D_os - D_ol.to(u.kpc) # Distance between lens and source

# Einstein radius (gravitational lensing angle) in milliarcseconds
theta_E_rad = np.sqrt((4*cst.G*M_lens)/(cst.c**2) * (D_ls/(D_ol*D_os))) * u.rad
theta_E = theta_E_rad.to(u.mas)
print(f"Einstein radius: {theta_E:.1f}")


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
    sep = np.linalg.norm(theta, axis=0).to(u.rad).value + 1e-12 # Angular separation

    return theta*(1 - (4*cst.G*mass)/(cst.c**2*D_L * sep**2) * (D_LS/D_L)) # Lensing formula

def get_axes(im):
    """Simplified axis generation - fixed scale in mas"""
    fieldOfView = 200000  # Field of view in mas # Equal to 3 arcminutes
    xs = np.linspace(-fieldOfView/2, fieldOfView/2, im.shape[0])
    ys = np.linspace(-fieldOfView/2, fieldOfView/2, im.shape[1])

    # xs = np.linspace(-10, 10, im.shape[0])  # Fixed 20mas FOV
    # ys = np.linspace(-10, 10, im.shape[1])  # Matches image aspect ratio
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
    #original = interpolate(interpolators, xv.value, yv.value)
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
    #xs, ys = get_axes(im, margin_factor = 8)
    xs, ys = get_axes(im, margin_factor = 8)
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

def animateWarping(im, outputPath = "microlensAnimation.gif", nFrames = 100):
    xs, ys = get_axes(im)
    interpolators = get_interpolators(im, xs, ys, interpolate_value=0.)
    extent = np.array([xs.min(), xs.max(), ys.min(), ys.max()])
    
    fig, ax = plt.subplots(figsize=(5, 5*im.shape[1]/im.shape[0]), dpi=512)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.axis('off')

    # mesh grid
    xv, yv = np.meshgrid(xs, ys) * u.mas
    imDisplay = ax.imshow(np.zeros_like(im), extent=extent, origin='lower')

    def moveFrame(frame):
        # fractional progress of animation
        fractional = frame / (nFrames - 1)

        # moves warping based on fractional progress
        xPos = extent[0] + fractional * (extent[1] - extent[0])
        yPos = 0  # Keep y fixed

        # Position array for black hole
        BH_centre = np.array([xPos, yPos]) * u.mas

        warped = apply_lensing(BH_centre, interpolators, xv, yv)
        imDisplay.set_data(warped)
        return [imDisplay]
    
    # calls moveFrame() nFrame times
    animatedLensing = animation.FuncAnimation(fig, moveFrame, frames=nFrames, blit=True)

    # animatedLensing defaults to using pillow so output file is a gif and not mp4
    animatedLensing.save(outputPath, fps=30)
    plt.close(fig)
#     print(f"Animation saved to: {outputPath}")

def animateBHLensing(im, iPosition, fPosition, outputPath='blackHoleSim.gif', nFrames=100):
    # Physical setup (adjusted for visible effect)
    # M_lens = 1000*u.M_sun  # More massive lens for visible effect
    # D_os = 8*u.kpc         # More reasonable distance
    # D_ol = 4*u.kpc
    
    xs, ys = get_axes(im)
    interpolators = get_interpolators(im, xs, ys, interpolate_value=0.)
    extent = np.array([xs.min(), xs.max(), ys.min(), ys.max()])

    fig, ax = plt.subplots(figsize=(10, 10*im.shape[1]/im.shape[0]), dpi=150)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.axis('off')

    xv, yv = np.meshgrid(xs, ys) * u.mas
    imDisplay = ax.imshow(im, extent=extent, origin='lower')

    # Red dot marker
    marker = plt.Circle((0, 0), radius=0.1, color='red', alpha=0.8, zorder=100)
    ax.add_patch(marker)

    # Convert positions to display coordinates
    init_ra, init_dec = iPosition[0][0].to(u.mas).value, iPosition[0][1].to(u.mas).value
    final_ra, final_dec = fPosition[0][0].to(u.mas).value, fPosition[0][1].to(u.mas).value

    # Scale positions to fit within image bounds
    all_ra = np.array([init_ra, final_ra])
    all_dec = np.array([init_dec, final_dec])
    
    scale_factor = min(
        (extent[1] - extent[0]) / (all_ra.max() - all_ra.min()) * 0.8,
        (extent[3] - extent[2]) / (all_dec.max() - all_dec.min()) * 0.8
    )
    
    def moveFrame(frame):
        fractional = frame / (nFrames - 1)
        
        # Calculate true BH position
        ra = init_ra + fractional * (final_ra - init_ra)
        dec = init_dec + fractional * (final_dec - init_dec)
        
        # Scale position for display
        ra_display = (ra - all_ra.mean()) * scale_factor
        dec_display = (dec - all_dec.mean()) * scale_factor
        
        # Update marker position
        marker.center = (ra_display, dec_display)
        
        # Apply lensing at display position (not TRUE position)
        BH_centre = np.array([ra_display, dec_display]) * u.mas
        warped = apply_lensing(BH_centre, interpolators, xv, yv)
        imDisplay.set_data(warped)

        return [imDisplay, marker]

    # Save animation
    writer = PillowWriter(fps=30)
    anim = animation.FuncAnimation(fig, moveFrame, frames=nFrames, blit=True)
    anim.save(outputPath, writer=writer)
    plt.close(fig)

# if __name__ == '__main__':
#     assert len(sys.argv) == 3, "Usage: python microlens_image.py <image> <output>"
    
#     # Load image from command line argument
#     im = load_image(sys.argv[1])
#     iPosition, fPosition = getProperMotion(n=30, years=1)

#     # warp_image(im, np.array([0, 0])*u.mas, sys.argv[2])
#     #animateWarping(im, outputPath=sys.argv[2], nFrames=100)
#     animateBHLensing(im, [iPosition[0]], [fPosition[0]], outputPath=sys.argv[2], nFrames=100)

if __name__ == '__main__':
    assert len(sys.argv) == 3, "Usage: python microlens_image.py <image> <output_folder>"
    print(theta_E)

    outputFolder = sys.argv[2]
    output_dir = Path(outputFolder)
    output_dir.mkdir(exist_ok=True)

    # Load image from command line argument
    im = load_image(sys.argv[1])
    print(f"Image loaded - shape: {im.shape}, dtype: {im.dtype}, range: [{im.min()}, {im.max()}]")
    plt.imshow(im)
    plt.title("Original Image Check")
    plt.savefig("image_check.png")
    print("Saved image_check.png - please verify this looks correct")
    
    # Get positions for all black holes
    iPosition, fPosition = getProperMotion(n=300, years=5)
    
    # Filter for BHs with |μ| ≥ 0.01 arcsec/yr
    highMotionBH = []
    for i in range(len(iPosition)):
        # Calculate proper motion from position change
        raD = (fPosition[i][0] - iPosition[i][0]).to(u.mas)
        decD = (fPosition[i][1] - iPosition[i][1]).to(u.mas)
        pm = np.sqrt(raD**2 + decD**2)  # Total proper motion in mas/yr
        if abs(pm.value) >= 10:  # 10 mas/yr = 0.01 arcsec/yr
            highMotionBH.append(i)

    # Generate GIFs for first 10 qualifying BHs
    for i, index in enumerate(highMotionBH[:100]):
        print(f"\nGenerating animation for BH {i+1} (index {index})")
        
        # Convert positions to Quantity objects in mas
        init_pos = (iPosition[index][0].to(u.mas), iPosition[index][1].to(u.mas))
        final_pos = (fPosition[index][0].to(u.mas), fPosition[index][1].to(u.mas))
        
        output_path = f"lensing{i+1}.gif"  # Names: lensing1.gif, lensing2.gif, etc.
        output_path = output_dir / f"lensing_{i+1}.gif"
        animateBHLensing(im, [init_pos], [final_pos], outputPath=str(output_path), nFrames=100)

    print(f"\nAll animations saved to: {output_dir.resolve()}")