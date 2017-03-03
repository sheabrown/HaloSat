# ============================================================================
# Functions for computing the observing window for the HaloSat X-ray satilite, 
# as well as misc. functions for tiling the sky with the 10.1 diameter beam.
# Used in the Astrophysical Machine Learning course at the University of Iowa
# https://astrophysicalmachinelearning.wordpress.com/ taught by Shea Brown
# Written by Shea Brown, shea-brown@uiowa.edu, and ... 
# https://sheabrownastro.wordpress.com/
# ============================================================================
import pyfits
import numpy as np
import healpy as hp
import numpy.ma as ma
import ephem
pi=np.pi
nan=np.nan


def getObsWindow(date,obs_window):
	# ===================================
	# date is a string with the 
	# Universal time, e.g. '2019/3/9 5:13'
	# obs_window is the radius(in degrees)
	# of the maximum deviation from the 
	# Sun's antipodal point that HaloSat
	# is allowed to observe. 
	# ===================================

	# Define the Sun object
	# ----------------------
	sun=ephem.Sun()

	# Define the observer
	# -------------------
	halosat = ephem.Observer()
	halosat.date = date

	# Get the location of the Sun
	# ----------------------------
	sun.compute(halosat)
	hours=ephem.hours
	deg=ephem.degrees
	
	# Print out the antipodal point 
	# -------------------------------
	print('On the data',halosat.date)
	print('The RA of the Sun is',sun.ra)
	print('The Dec of the Sun is',sun.dec)
	print('So the HaloSat will observe')
	win_ra  = hours(sun.ra+pi).norm
	win_dec = deg(sun.dec*(-1.))
	print('      RA ',win_ra) 
	print('      DEC',win_dec)

	# Read map of the galaxy for plotting the window
	# ----------------------------------------------
	map = hp.read_map('ROSAT_Total_hp_filt.fits')
	cone_radius=obs_window/57.2958 # convert to radians
	
	# Convert to Galactic coordinates
	# -------------------------------
	ap = ephem.Equatorial(win_ra, win_dec, epoch='2000')
	sp = ephem.Equatorial(sun.ra, sun.dec, epoch='2000')
	g = ephem.Galactic(ap)
	s = ephem.Galactic(sp)
	print(g.lon,g.lat)
	vec=hp.ang2vec(g.lon*57.4,g.lat*57.2958,lonlat=True)
        local=hp.query_disc(512,vec,cone_radius)
	map[local]=nan
	hp.mollview(map,max=6000,title='Date: '+str(halosat.date)+', Window Center: [RA '+str(win_ra)+', DEC ' + str(win_dec)+']')
	hp.projtext(g.lon*57.2958+20, g.lat*57.2958, 'Obs Window', lonlat=True, coord='G',fontsize=20)
	hp.projscatter(s.lon*57.2958, s.lat*57.2958, lonlat=True, s=55,coord='G',color='yellow')
	hp.projtext(s.lon*57.2958-5, s.lat*57.2958, 'Sun', lonlat=True, coord='G',fontsize=20)

def getRandomPointings(nside,N):
	npix=hp.nside2npix(nside)
	return np.random.randint(npix,size=N)

def plotPointings(nside,pointings):
	npix=hp.nside2npix(nside)
	map=np.zeros(npix)
	rad=5.05
	rad=rad/57.3 #convert to radians
	for x in range(0,len(pointings),1):
                pix=hp.pix2vec(nside,pointings[x])
                local=hp.query_disc(nside,pix,rad)
                map[local]=1.0
	hp.mollview(map)


# Try it out
# ------------------
#pointings=getRandomPointings(128,300)
#print(len(pointings))
#plotPointings(128,pointings)
