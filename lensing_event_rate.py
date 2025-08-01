import csv
import math
import numpy as np
import astropy.units as u
from astropy.time import Time
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, CartesianRepresentation, Galactocentric, ICRS, Angle
import os

def getProperMotion(n = 100000, years = 5):
    iPosition = []
    fPosition = []
    bhCount = 0
    selected = 0
    distances = []
    proper_motions = []

    with open(r'C:\Users\User\Downloads\kicked_remnants_igoshev_young_7.8_DC_integrated_final.csv', 'r') as file:
        lines = csv.reader(file)
        next(lines)

        count = 0
        for line in lines:
            if line[13].strip() == "Black Hole":

                bhCount += 1
                # Position in kpc and velocity in km/s
                xkPC, ykPC, zkPC = float(line[2]), float(line[3]), float(line[4])
                vx, vy, vz = float(line[5]), float(line[6]), float(line[7])

                initial_pos = SkyCoord(
                    x=xkPC * u.kpc,
                    y=ykPC * u.kpc,
                    z=zkPC * u.kpc,
                    frame='galactocentric'
                )

                posICRS = initial_pos.transform_to('icrs')
                dist = posICRS.distance.to(u.pc).value # in pc
                ra = posICRS.ra.to(u.deg).value
                dec = posICRS.dec.to(u.deg).value

                # Position array in parsecs
                # pos = np.array([x.to(u.pc).value, y.to(u.pc).value, z.to(u.pc).value])
                #pos = np.array([x.value, y.value, z.value])
                #pos = np.array([x, y, z])
                #vel = np.array([vx, vy, vz])

                vRadial = float(line[17])
                vTransverse = float(line[19])
                vTotal = float(line[20])

                v_total_calculated = np.sqrt(vx**2 + vy**2 + vz**2)
                v_trans_calculated = np.sqrt(v_total_calculated**2 - vRadial**2)


                difference = v_trans_calculated - vTransverse


                # Proper mu motion in arcsec/year
                mu = (vTransverse / (4.74 * dist)) * u.arcsec/u.yr # 4.74 constant from converting kpc to km and mas to rad/s
                #print(f"Proper Motion: {mu}") if mu > 0.01 * u.arcsec/u.yr else print(None)
                #mu = vTransverse / (4.74 * dist.to(u.kpc).value)
                #if mu < 1 * u.arcsec/u.yr:  # Skip nearly static BHs
                   #continue                

                # Position in cartesian coordinates
                cart = CartesianRepresentation(x=xkPC * u.kpc, y=ykPC * u.kpc, z=zkPC * u.kpc)
                coord_cart = SkyCoord(cart,
                                    frame='galactic',
                                    representation_type='cartesian',
                                    obstime=Time('2000-01-01'))
                
                # Convert to spherical galactic coordinates (Refer to astropy API)
                l = coord_cart.spherical.lon
                b = coord_cart.spherical.lat
                distance = coord_cart.spherical.distance

                pm_l_cosb = (mu / np.sqrt(2))
                pm_b = (mu / np.sqrt(2))

                #SkyCoord object with proper motion
                coord = SkyCoord(l=l, b=b, distance=distance,
                                pm_l_cosb=pm_l_cosb, pm_b=pm_b,
                                frame='galactic',
                                obstime=Time('2000-01-01'))
                
                
                age = float(line[8])
                # Evolve across sky for 1 year
                coordEvolved = coord.apply_space_motion(new_obstime=Time('2001-01-01') + age * u.yr + years * u.yr)


                # Convert SkyCoord to (RA, Dec) in mas
                if mu >= 0.05 * u.arcsec/u.yr or mu <= -0.05 * u.arcsec/u.yr:
                    iPosition.append((coord.icrs.ra, coord.icrs.dec))  # Astropy objects
                    fPosition.append((coordEvolved.icrs.ra, coordEvolved.icrs.dec))
                    distances.append(dist / 1000)  # Convert to kpc for plotting
                    proper_motions.append(mu.value)  # Extract numerical value

                count += 1
                if count >= n:
                    break
                # Print positions in RA/Dec
                #print("Initial Position:", coord.icrs.to_string('hmsdms'))
                #print("Position After 1 year:", coordEvolved.icrs.to_string('hmsdms'))
                    
                if mu >= 0.05 * u.arcsec/u.yr or mu <= -0.05 * u.arcsec/u.yr:
                     selected += 1
                     print(f"\nBlack Hole {selected} (μ = {mu})")
                     print(f"RA: {int(coord.icrs.ra.hms.h)}h {int(coord.icrs.ra.hms.m)}m {coord.icrs.ra.hms.s}s")
                     print(f"Dec: {int(coord.icrs.dec.degree)}° {int(coord.icrs.dec.arcminute)}' {coord.icrs.dec.arcsecond}\"")
                     print(f"Final RA: {int(coordEvolved.icrs.ra.hms.h)}h {int(coordEvolved.icrs.ra.hms.m)}m {coordEvolved.icrs.ra.hms.s}s")
                     print(f"Final Dec: {int(coordEvolved.icrs.dec.degree)}° {int(coordEvolved.icrs.dec.arcminute)}' {coordEvolved.icrs.dec.arcsecond}\"")
                    

        #print(f"\nFound {selected} BHs with μ ≥ 0.05 arcsec/yr (out of {bhCount} total)")
        print(f"\nFound {selected} BHs with distance < 1 kpc (out of {bhCount} total)")
        plotHistograms(distances, proper_motions)
        return iPosition, fPosition
    

def plotHistograms(distances, proper_motions):
    
    # Proper motion to absolute value
    proper_motions_abs = [abs(pm) for pm in proper_motions]
    
    # Figures with subplts
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # Black hole Distance distribution
    ax1.hist(distances, bins=20, alpha=0.7, color='blue', edgecolor='black')
    ax1.set_xlabel('Distance from Earth (kpc)')
    ax1.set_ylabel('Number of Black Holes')
    ax1.set_title('Distance Distribution of Black Holes with Proper Motion |μ| ≥ 0.05')
    ax1.grid(True, alpha=0.3)
    
    # Proper motion distribution
    ax2.hist(proper_motions_abs, bins=20, alpha=0.7, color='red', edgecolor='black')
    ax2.set_xlabel('Proper Motion (arcsec/yr)')
    ax2.set_ylabel('Number of Black Holes')
    ax2.set_title('Distribution of Black Holes with Proper Motion |μ| ≥ 0.05')
    ax2.set_yscale('log')  # Log scale to take into account large y range
    ax2.grid(True, alpha=0.3)
    
    # 2D histogram of black hole distance distribution
    ax3.hist2d(distances, proper_motions_abs, bins=50, cmap='viridis')
    ax3.set_xlabel('Distance from Earth (kpc)')
    ax3.set_ylabel('Proper Motion (arcsec/yr)')
    ax3.set_title('Distance vs Proper Motion |μ| ≥ 0.05')
    ax3.set_yscale('log')
    
    # Scatter plot
    ax4.scatter(distances, proper_motions_abs, alpha=0.6, s=20, color='purple')
    ax4.set_xlabel('Distance from Earth (kpc)')
    ax4.set_ylabel('Proper Motion (arcsec/yr)')
    ax4.set_title('Distance vs Proper Motion |μ| ≥ 0.05')
    ax4.set_yscale('log')
    ax4.set_ylim(bottom = None, top = 1)
    ax4.grid(True, alpha=0.3)

    
    plt.tight_layout()
    plt.show()
    
    # Some stats
    print(f"\nStatistics for Black Holes with |μ| ≥ 0.1 arcsec/yr:")
    print(f"  Number of high proper motion BHs: {len(distances)}")
    print(f"\nDistance Statistics:")
    print(f"  Mean distance: {np.mean(distances):.2f} kpc")
    print(f"  Median distance: {np.median(distances):.2f} kpc")
    print(f"  Min distance: {np.min(distances):.2f} kpc")
    print(f"  Max distance: {np.max(distances):.2f} kpc")
    
    print(f"\nProper Motion Statistics:")
    print(f"  Mean proper motion: {np.mean(proper_motions_abs):.4f} arcsec/yr")
    print(f"  Median proper motion: {np.median(proper_motions_abs):.4f} arcsec/yr")
    print(f"  Min proper motion: {np.min(proper_motions_abs):.4f} arcsec/yr")
    print(f"  Max proper motion: {np.max(proper_motions_abs):.4f} arcsec/yr")

def getBHDistances(n = 100000, years = 5):
    bhPositions = []    

    with open(r'C:\Users\User\Downloads\kicked_remnants_igoshev_young_7.8_DC_integrated_final.csv', 'r') as file:
        lines = csv.reader(file)
        next(lines)

        count = 0
        for line in lines:
            if line[13].strip() == "Black Hole":

                # Position in kpc and velocity in km/s
                xkPC, ykPC, zkPC = float(line[2]), float(line[3]), float(line[4])
                vx, vy, vz = float(line[5]), float(line[6]), float(line[7])
                #print(xkPC, ykPC, zkPC)
                #print(type(xkPC))

                bhPositions.append((xkPC, ykPC, zkPC))

                count += 1
                if count > n:
                    break

    radius = 8.2 * u.kpc
    z = 0 * u.kpc
    nEarths = 1000

    # Array of angles/positions around the galaxy centre
    angles = np.linspace(0, 2 * np.pi, nEarths)

    output_dir = "histograms"
    os.makedirs(output_dir, exist_ok=True)

    bhCoords = [
    SkyCoord(x=x * u.kpc, y=y * u.kpc, z=z * u.kpc, frame="galactocentric")
    for x, y, z in bhPositions
]


    for i, theta in enumerate(angles): 
        bhDistances = []
        earthCoord = SkyCoord(
        x=radius * np.cos(theta),
        y=radius * np.sin(theta),
        z=z,
        frame="galactocentric"
    )
        bhDistances = [earthCoord.separation_3d(bh).to(u.kpc).value for bh in bhCoords]
        bhDistances = [d for d in bhDistances if d <= 1.0]

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(bhDistances, bins=25)
        ax.set_title(f"Black Hole Distances for simulated Earth {i} at angle {theta}")
        ax.set_xlabel("Distance from Earth (kpc)")
        ax.set_ylabel("Number of Black Holes")
        ax.set_xlim(0, 1)
        ticks = np.arange(0, 1.1, 0.1)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{tick:.1f}" for tick in ticks])
        ax.grid(True)
        
        filename = f"Simulated Earth {i:04}.png"
        filepath = os.path.join(output_dir, filename)
        fig.savefig(filepath)
        plt.close()

                
                    