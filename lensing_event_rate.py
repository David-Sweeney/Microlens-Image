import csv
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, CartesianRepresentation

def getProperMotion(n = 10):
    iPosition = []
    fPosition = []

    with open(r'C:\Users\User\Downloads\kicked_remnants_igoshev_young_7.8_DC_integrated_final.csv', 'r') as file:
        lines = csv.reader(file)
        next(lines)

        count = 0
        for line in lines:
            if line[13].strip() == "Black Hole":

                # Position in kpc and velocity in km/s
                xkPC, ykPC, zkPC = float(line[2]), float(line[3]), float(line[4])
                vx, vy, vz = float(line[5]), float(line[6]), float(line[7])

                # Position to astropy
                x = xkPC * u.kpc
                y = ykPC * u.kpc
                z = zkPC * u.kpc

                # Position array in parsecs
                pos = np.array([x.to(u.pc).value, y.to(u.pc).value, z.to(u.pc).value])
                vel = np.array([vx, vy, vz])

                vRadial = float(line[17])
                vTransverse = float(line[19])
                vTotal = float(line[20])

                # Distance with normalised position vector
                dist = np.linalg.norm(pos)  # in pc

                # Proper mu motion in arcsec/year
                mu = vTransverse / (4.74 * dist) # 4.74 from converting kpc to km and mas to rad/s

                age = float(line[8])
                

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

                pm_l_cosb = (mu / np.sqrt(2)) * u.arcsec/u.yr
                pm_b = (mu / np.sqrt(2)) * u.arcsec/u.yr
                #pm_l_cosb = round(mu * u.arcsec/u.yr / np.sqrt(2), 2)
                #pm_b = round(mu * u.arcsec/u.yr / np.sqrt(2), 2)

                #SkyCoord object with proper motion
                coord = SkyCoord(l=l, b=b, distance=distance,
                                pm_l_cosb=pm_l_cosb, pm_b=pm_b,
                                frame='galactic',
                                obstime=Time('2000-01-01'))

                # Evolve across sky for 1 year
                coord1yL = coord.apply_space_motion(new_obstime=Time('2001-01-01') + age * u.yr + 1 * u.yr)

                iPosition.append(coord.icrs)
                fPosition.append(coord1yL.icrs)
                count += 1

                if count >= n:
                    break

            return iPosition, fPosition
                
                    # Print positions in RA/Dec
                    # print("Initial Position:", coord.icrs.to_string('decimal'))
                    # print("Position After 1 year:", coord1yL.icrs.to_string('decimal'))