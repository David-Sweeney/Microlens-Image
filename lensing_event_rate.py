import csv
import math
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, CartesianRepresentation

def getProperMotion(n = 300, years = 5):
    iPosition = []
    fPosition = []
    bhCount = 0
    selected = 0

    with open(r'C:\Users\User\Downloads\kicked_remnants_igoshev_young_7.8_DC_integrated_final.csv', 'r') as file:
        lines = csv.reader(file)
        next(lines)

        count = 0
        for line in lines:
            if line[13].strip() == "Black Hole":

                bhCount += 1
                # Position in kpc and velocity in km/s
                xkPC, ykPC, zkPC = float(line[2]), float(line[3]), float(line[4])
                #vx, vy, vz = float(line[5]), float(line[6]), float(line[7])

                # Position to astropy
                x = xkPC * u.kpc
                y = ykPC * u.kpc
                z = zkPC * u.kpc

                # Position array in parsecs
                pos = np.array([x.to(u.pc).value, y.to(u.pc).value, z.to(u.pc).value])
                #pos = np.array([x, y, z])
                #vel = np.array([vx, vy, vz])

                vRadial = float(line[17])
                vTransverse = float(line[19])
                vTotal = float(line[20])

                # Distance with normalised position vector
                dist = np.linalg.norm(pos)  # in pc

                # Proper mu motion in arcsec/year
                #mu = vTransverse / (4.74 * dist) 
                mu = (vTransverse / (4.74 * dist)) * u.arcsec/u.yr # 4.74 from converting kpc to km and mas to rad/s
                #print(f"Proper Motion: {mu}") if mu > 0.01 * u.arcsec/u.yr else print(None)
                #mu = vTransverse / (4.74 * dist.to(u.kpc).value)
                #if mu < 1 * u.arcsec/u.yr:  # Skip nearly static BHs
                   #continue                

                # Position in cartesian coordinates
                cart = CartesianRepresentation(x=x, y=y, z=z)
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
                #pm_l_cosb = round(mu * u.arcsec/u.yr / np.sqrt(2), 2)
                #pm_b = round(mu * u.arcsec/u.yr / np.sqrt(2), 2)

                #SkyCoord object with proper motion
                coord = SkyCoord(l=l, b=b, distance=distance,
                                pm_l_cosb=pm_l_cosb, pm_b=pm_b,
                                frame='galactic',
                                obstime=Time('2000-01-01'))
                
                age = float(line[8])
                # Evolve across sky for 1 year
                coordEvolved = coord.apply_space_motion(new_obstime=Time('2001-01-01') + age * u.yr + years * u.yr)

                # iPosition.append(coord.icrs)
                # fPosition.append(coord1yL.icrs)

                # Convert SkyCoord to (RA, Dec) in mas
                iPosition.append((coord.icrs.ra, coord.icrs.dec))  # Astropy objects
                fPosition.append((coordEvolved.icrs.ra, coordEvolved.icrs.dec))
                #iPosition.append((coord.icrs.ra.degree, coord.icrs.dec.degree)) # Regular float values
                #fPosition.append((coordEvolved.icrs.ra.degree, coordEvolved.icrs.dec.degree))
                #iPosition.append((coord.icrs.ra.to(u.mas).value, coord.icrs.dec.to(u.mas).value))
                #fPosition.append((coordEvolved.icrs.ra.to(u.mas).value, coordEvolved.icrs.dec.to(u.mas).value))

                count += 1
                if count >= n:
                    break
                # Print positions in RA/Dec
                #print("Initial Position:", coord.icrs.to_string('hmsdms'))
                #print("Position After 1 year:", coordEvolved.icrs.to_string('hmsdms'))
                    
                if mu >= 0.01 * u.arcsec/u.yr or mu <= -0.01 * u.arcsec/u.yr:
                     selected += 1
                     print(f"\nBlack Hole {selected} (μ = {mu})")
                     print(f"RA: {int(coord.icrs.ra.hms.h)}h {int(coord.icrs.ra.hms.m)}m {coord.icrs.ra.hms.s}s")
                     print(f"Dec: {int(coord.icrs.dec.degree)}° {int(coord.icrs.dec.arcminute)}' {coord.icrs.dec.arcsecond}\"")
                     print(f"Final RA: {int(coordEvolved.icrs.ra.hms.h)}h {int(coordEvolved.icrs.ra.hms.m)}m {coordEvolved.icrs.ra.hms.s}s")
                     print(f"Final Dec: {int(coordEvolved.icrs.dec.degree)}° {int(coordEvolved.icrs.dec.arcminute)}' {coordEvolved.icrs.dec.arcsecond}\"")
                    
        
        #print(iPosition)
        #print(fPosition)
        print(f"\nFound {selected} BHs with μ ≥ 0.01 arcsec/yr (out of {bhCount} total)")
        return iPosition, fPosition
    
#getProperMotion(n=10, years=1)
                
                    