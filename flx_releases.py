#!/usr/bin/env python
# flx_releases.py
from numpy import *                             # Numerics
import asciitable                               # ASCII IO
from datetime import *                          # Manipulate dates
import warnings
import re                                       # regex tools
import math                                     # Math tools

# Global variables:
pathOptions = "./options"   # path to COMMAND and RELEASES files
speciesIndex=15             # table entry
summitHeight = 1.55         # km asl
centerLat   = 63.63         # Latitude of source (center)
centerLon   =-19.63         # Longitude of source (center)
npart       = 100000        # Number of particles released per hour
mass_fraction = 1.0     # For multiple spec. runs; fraction of spec./total mass
mixing_ratio= 1.0e+2        # Desired mixing ratio in the eruption column[kg/kg]
                            # This is not the actual mixing ratio in the plume
                            # as it is way too high, to be buoyant.
                            # Setting it to a "realistic" value will cause a 
                            # plume spread of ~1000 km...
                            # However, not all ash in the column is suspended,
                            # some will simply be there due to inertial upward
                            # motion.

# Notes on global variables:
# mass_fraction     - fraction of total mass allocated to the species simulated
#                     in the current run. This should probably be 1 for single 
#                     species runs.
# conc_at_sea_level - We calculate the concentration throughout the plume
#                     based on a virtual concentration at sea level. What this
#                     should be is a bit unclear, sorry.  We should make some
#                     tests on this, it should probably depend on total mass or
#                     eruption column height in some way. We should try some
#                     different value and compare it to estimates of ash load at
#                     the top of eruption columns.

#TODO
# def setTime(simstart):

class Sources:
    def __init__(self,length=1,basename="default"):
        if type(length) != int:
            raise Exception("Wrong argument, expected",int,\
                            "but got",type(length))
        self.n  = length
        self.t0 = [datetime(1900,1,1) for i in range(length)]   # start time
        self.t1 = [datetime(1900,1,1) for i in range(length)]   # end time
        self.x0 = zeros(length)     # Longitude lower boundary
        self.x1 = zeros(length)     # Longitude upper boundary
        self.y0 = zeros(length)     # Latitude lower boundary
        self.y1 = zeros(length)     # Latitude upper boundary
        self.z0 = zeros(length)     # Altitude lower boundary
        self.z1 = zeros(length)     # Altitude upper boundary
        self.aglasl = 2             # Altitude in m asl (1=agl, not sup.)
        self.mass   = zeros(length) # Total mass emitted (kg)
        self.npart  = zeros(length) # Number of model particles emitted
        self.name   = [''.join([basename,"_z",str(i)]) for i in range(length)]

    def metres2degrees(self,metres):
        # Get distance in latitude,longitude from given distance in m.
        # This method assumes a spherical Earth and that the distance 
        # from centerLat is small. The method should be accurate to about
        # 1 % for latitude, and >1 % for Longitude (near 1 % at equator
        # but increases towards the poles, about 5 % at N63-N64)
        # Also check out: http://www.csgnetwork.com/degreelenllavcalc.html
        rEarth  = 6.371e6           # mean radius
        oEarth  = 2*math.pi*rEarth  # circumference
        m_per_deg_lat = oEarth/360
        m_per_deg_lon = m_per_deg_lat*math.cos(math.radians(centerLat))
        dLat    = metres/m_per_deg_lat
        dLon    = metres/m_per_deg_lon

        return dLon,dLat
        
    def std_atm(self,alt):
        # Function should accept 1-d arrays or scalars
        # Takes altitude in km as input
        # Returns density, pressure and temperature ratios (in that order)
        # To get the absolute of any value, multiply it with the ground level
        # value, e.g. T(z) = ratio(z)*T0

        # CONSTANTS
        rearth  = 6369.0        # Radius of the Earth
        gmr     = 34.163195     # Hydrostatic constant
        ntab    = 8             # No. of vertical section in std atm
        #
        
        # STANDARD ATMOSPHERE TABLES
        htab = array([0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852],\
                     dtype=float)
        ttab = array([288.15, 216.65, 216.65, 228.65, \
                      270.65, 270.65, 214.65, 186.946],\
                     dtype=float)
        ptab = array([1.0000000E-0, 2.2336110E-1, 5.4032950E-2, 8.5666784E-3,\
                     1.0945601E-3, 6.6063531E-4, 3.9046834E-5, 3.68501E-6],\
                     dtype=float)
        gtab = array([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0],dtype=float)
        # htab - geopotential height
        # ttab - temperature
        # ptab - pressure
        # gtab - lapserate (g = gamma) 
        #

        # BEGIN
        alt*rearth
        height = alt*rearth/(alt+rearth) # convert geometric to geopotential height

        theta = zeros(shape(height)[:])     # predefine some arrays
        delta = zeros(shape(height)[:])
        for i,h in ndenumerate(height):
            ilo = max(nonzero(htab<=h))[0] # which table entry to use
            tgrad = gtab[ilo]           # temperature gradient in layer
            tbase = ttab[ilo]           # temperature at bottom of layer
            deltah  = h-htab[ilo]       # distance from bottom of layer
            tlocal  = tbase+tgrad*deltah    # temperature at given height
            theta[i] = tlocal/ttab[0]       # temperature ratio (local/sea level)
            if tgrad == 0:
                delta[i] = ptab[ilo]*exp(-gmr*deltah/tbase) # pressure ratio
            else:
                delta[i] = ptab[ilo]*(tbase/tlocal)**(gmr/tgrad) # pressure ratio

        sigma = delta/theta     # Density ratio

        return sigma,delta,theta

    def calcLayers(self,hgts,mass,scaleFactors):
        #TODO: check arguments
        # input hgts in [m]
        # sets up the source layers based on total column height
        # hgts is a list of vertical layer tops
        # p is an array with scale factors (from probability density function)
        # p will be normalized so that sum(p)=1
        # We also require the global variable summitHeight to determine the
        # bottom of the lowermost layer.
        # Standard atmospheric data will be used to determine horizontal
        # spread of the plume
        
        # We want mixing ratio to be about constant throughout the eruption
        # column, this is the same as a constant ratio between the atmospheric
        # density and concentration of ash should be constant throughout the
        # column:
        #   [ash conc]/[air density] = [MassMix] = constant
        # If we know the total mass emissions rate [mer] at a thin vertical
        # segment, we can express the volume of that segment as follows:
        #   dV = X(z)*Y(z) dz = [residence time]*[mer]/[MassMix]/[air density]
        # Where the residence time is the average time it takes for newly
        # released particle to leave the source volume laterally.
        # The residence time should be rev. proport. to the horiz. wind speed
        # [residence time] = 1/U * X(z)/2
        # Thus, if we know the air density and the mass mixing ratio we
        # can determine a theoretical value for the horizontal spread of
        # the plume: 
        #   X(z) = 1/(2*[MassMix]) * [mer]/([density] U dz] (1)
        #
        # A limitation of this method is that we do not account the mass 
        # loss due to fallout of larger clasts in the lower portion of the 
        # eruption column, however, as we are mainly interested in the
        # long range dispersion, we probably don't even want to consider
        # such large clasts.  This assumption gets less valid for larger
        # eruption columns(?).
        #
        # It would be preferable to assign a cylindrical shape of each
        # segment but the current version of flexpart cannot support this
        # so we have to make do with box-shaped plume segments

        self.z1[:]  = [h*1000 for h in hgts[:]]       # Top of layers
        self.z0[1:] = [h*1000 for h in hgts[:-1]]     # Bottom of layers
        self.z0[0]  = summitHeight*1000  # Bottom of lowest layer

        # Make an estimate of the air density in the ambient atmosphere
        rho_air = zeros(shape(hgts)[0])   # Desired ash concentration at heights
        u_ave = zeros(shape(hgts)[0])     # Wind speed at heights
        for ih in xrange(shape(hgts)[0]):
            # Create a sample of heights in the layer (to improve accuracy)
            # TODO: Number of samples should depend on layer span and/or
            #       second derivative of density with regard to z
            h = linspace(self.z1[ih],self.z0[i],10)/1000    # in km !!
            # Calculate standard atmospheric parameters for each sample
            sigmas,deltas,thetas = self.std_atm(h)
            # Here, sigma is the ratio between sea level density and density
            # the given altitude.  Thus it is also the ratio between ash mass
            # MR relative a virtual MR at sea level.
            # We will use the mean of the ratio to get the bulk air density
            rho_air[ih] = mean(sigmas)

            # We currently just set wind speed to a constant:
            u_ave[ih] = 10      # m/s

        # Normalize mass scaleFactors
        sfm = scaleFactors/sum(scaleFactors)        # for Mass

        # For each layer (slice) of the source
        print 'z, dx in metres:'   # output from loop, mainly for verifying result
        for ih in xrange(shape(hgts)[0]):

            # The mass emitted at a given layer is the total emission by the
            # volcano over the given time period time the fraction allocated to
            # the current species times the scale factor for the current level.
            self.mass[ih] = mass*mass_fraction*sfm[ih] 

            # The duration of the eruption phase is used to determine number of
            # particles released the shape of the source cloud
            diff = self.t1[ih]-self.t0[ih]      # duration as timedelta and...
            dt_hours = float(diff.days*24)+float(diff.seconds)/3600.0 #number
            dt_seconds = dt_hours*3600.0

            #print dt,diff.seconds
            if dt_hours <= 0.0:
                print "t0,t1:",self.t1[ih],"\t",self.t0[ih]
                print "dt_hours:",dt_hours
                raise Exception("Duration of phase is too short, make sure"\
                                "t0 and t1 are properly defined")

            # Calculate horizontal spread according to (1)
            # this should depend on the sea level virtual concentration (or some
            # predetermined mixing ratio?), air density at the given layer and
            # the mass emitted at the given layer as well as the duration of the
            # current phase.
            dz = self.z1[ih]-self.z0[ih]        #depth of layer
            dx_metres = 0.5/mixing_ratio*\
                    self.mass[ih]/(dt_seconds*rho_air[ih]*u_ave[ih]*dz)
            print self.z1[ih], dx_metres
            # Convert horizontal spread to degrees lat,lon
            dx,dy = self.metres2degrees(dx_metres)
            self.x0[ih] = centerLon-dx/2
            self.x1[ih] = centerLon+dx/2
            self.y0[ih] = centerLat-dy/2
            self.y1[ih] = centerLat+dy/2


            self.npart[ih] = npart/self.n*dt_hours # N particles to be released
            

    def fPrint(self,fid):
        # Appends source information to given file
        # All sources in the class will be written to the file
        if type(fid) != file:
            raise Exception("Wrong argument, expected",file,\
                            "but got",type(fid))

        for i in range(self.n):
            fid.write(datetime.strftime(self.t0[i],"%Y%m%d %H%M%S\n"))
            fid.write("________ ______            "\
                      "i8,1x,i6 Beginning date and time of release\n\n")
            fid.write(datetime.strftime(self.t1[i],"%Y%m%d %H%M%S\n"))
            fid.write("________ ______            "\
                      "i8,1x,i6 Ending date and time of release\n\n")
            fid.write("%9.4f\n" %self.x0[i])
            fid.write("____.____                  "\
                      "f9.4  Longitude [DEG] of lower left corner\n\n")
            fid.write("%9.4f\n" %self.y0[i])
            fid.write("____.____                  "\
                      "f9.4  Latitude [DEG] of lower left corner\n\n")
            fid.write("%9.4f\n" %self.x1[i])
            fid.write("____.____                  "\
                      "f9.4  Longitude [DEG] of upper right corner\n\n")
            fid.write("%9.4f\n" %self.y1[i])
            fid.write("____.____                  "\
                      "f9.4  Latitude [DEG] of upper right corner\n\n")
            fid.write("%9d\n" %self.aglasl)
            fid.write("_________                  "\
                      "i9    1 for m above ground, 2 for m above sea level\n\n")
            fid.write("%9.3f\n" %self.z0[i])
            fid.write("_____.___                  "\
                      "f10.3 Lower z-level (in m agl or m asl)\n\n")
            fid.write("%9.3f\n" %self.z1[i])
            fid.write("_____.___                  "\
                      "f10.3 Upper z-level (in m agl or m asl)\n\n")
            fid.write("%9d\n" %self.npart[i])
            fid.write("_________                  "\
                      "i9    Total number of particles to be released\n\n")
            str = "%9.4E" %self.mass[i]         # Exponential output requires
            str = ''.join([str[:7],str[8:],"\n"]) # some tweaking to get rid of
            fid.write(str)                      # the "+"-sign after the "E"
            fid.write("_.____E__                  "\
                      "e9.4  Total mass emitted\n\n")
            fid.write("%s\n" %self.name[i])
            fid.write("________________________________________   "\
                      "character*40 comment\n")
            fid.write("+++++++++++++++++++++++++++++++++++++++++++++++++++++++"\
                      "++++++++++++++++++\n")
            

def writeReleaseHeader(fid):
    # Writes the file header
    if type(fid) != file:
        raise Exception("Wrong argument, expected",file,\
                        "but got",type(fid))
    h="*************************************************************************\n"
    v="*                                                                       *\n"
    fid.write(h)
    fid.write(v)
    fid.write(v)
    fid.write(v)
    fid.write("*    Input file for the Lagrangian particle dispersion model"\
              "FLEXPART    *\n")
    fid.write("*                        Please select your options"\
              "                     *\n")
    fid.write(v)
    fid.write("*                  FILE AUTO-GENERATED BY FLX_RELEASES.PY"\
              "               *\n")
    fid.write(v)
    fid.write(h)
    h="+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
    fid.write(h)
    fid.write("  1\n")
    fid.write("___                        i3    "\
              "Total number of species emitted\n\n")
    fid.write("%3d\n" %speciesIndex)
    fid.write("___                        i3    "\
              "Index of species in file SPECIES\n\n")
    h="=========================================================================\n"
    fid.write(h)

    

def readSimTime():
    # Read simulation start and end times from COMMAND
    fname   = '/'.join([pathOptions,'COMMAND']) # path to COMMAND file
    pattern = re.compile(r'\d{8} \d{6}')  # create pattern YYYYMMDD HHMMSS
    fid     = open(fname)                       # open file
    datestr = ['xxxxxxxx xxxxxx','xxxxxxxx xxxxxx'] # will store dates
    #print datestr
    i = 0
    print "Searching for simulation start and end in file:",fname
    for line in fid:
        match=pattern.search(line)              # Search for pattern
        if match:                               # Found match?
            print 'Match number',i,'is:'
            print match.group()
            datestr[i] = match.group()          # Store match in datestr
            i=i+1
            if i>2:                             # Too many matches in file?
                raise Exception('Too many dates found in COMMAND, terminating')
    
    # Convert to Python date
    date0_py  = datetime.strptime(datestr[0],"%Y%m%d %H%M%S")   # sim start
    date1_py  = datetime.strptime(datestr[1],"%Y%m%d %H%M%S")   # sim end
    return date0_py,date1_py


#TODO def readColHeight
# Check for error codes:
# Warn when encountering error:
#   9999: Scan missing
#   2000: Plume obscured by precipitation clouds
#   1000: Plume below detection limit
def readColHeight():
    minhgt= 2.5                         # Assumed minimum detection limit
    fname = 'column_heights.txt'        # File to read

    table = asciitable.read(fname)      # Read file
    date0 = table.col1                  # extract start date
    time0 = table.col2                  # extract start time
    date1 = table.col3                  # extract end date
    time1 = table.col4                  # extract end time
    hgt   = table.col5                  # extract height
    mass  = table.col6                  # extract mass release
    name  = table.col7                  # extract name/comment
    #TODO check if column 6 is missing, use power law to calc
    # mass if that is the case.
    for i in range(hgt.shape[0]):       # For each column height
    #for h in hgt:
        h=hgt[i]
        if h >=999:                     # Assume error
            msg="column height greater than 999 km found, assuming error code"
            warnings.warn(msg)
            print "Row: %i,     Height: %f",i,h
            raise Exception("NOT SUPPORTED YET, ABORTING")
    # Bundle date and time in same string
    datestr0 = [date0[i]+time0[i] for i in range(date0.shape[0])]
    datestr1 = [date1[i]+time1[i] for i in range(date1.shape[0])]
    # Convert string date/time to python date format
    date0_py  = [datetime.strptime(t,"%Y-%m-%d%H:%M:%S") for t in datestr0]
    date1_py  = [datetime.strptime(t,"%Y-%m-%d%H:%M:%S") for t in datestr1]
    # Make time relative to first start date:
    dateBase = date0_py[0]
    date0_py = [d-dateBase for d in date0_py]
    date1_py = [d-dateBase for d in date1_py]
    return date0_py,date1_py,hgt,mass,name
        

def poisson(nu,mu):
    # Returns a Poisson probability density function
    # from the expected mean value mu[0] at the given
    # counts nu[1]
    if type(mu) != int:
        raise Exception("Wrong argument, expected",int,\
                        "but got",type(mu))

    p = [ exp(-mu)*power(mu,i)/math.factorial(i) for i in nu ]
    return p


#TODO MAIN
# readSimTime
# readColHeight
# for h in colHeight:
#   s[i]=source(h)
date0,date1,hgt,mass,name=readColHeight()    # date0 & date1 are relative simstart
simstart,simend=readSimTime()
release_start   = [simstart+i for i in date0]
release_end     = [simstart+i for i in date1]

nLayers = 10                    # Number of plume slices (separate sources)
ilayer = range(nLayers)         # Layer numbers (0 -- nLayers-1)
p = flipud(poisson(ilayer,3))   # Prob. dens. func. used for alloc. emissions
                                # We invert this since we want an "upside-down"
                                # Poisson shape to represent the top-hat style
                                # eruption.

#fname   = '/'.join([pathOptions,'RELEASES']) # path to COMMAND file
#fout = open(fname,'w')
fout = open('/'.join([pathOptions,'RELEASES']),'w')      # Open output file
writeReleaseHeader(fout)    # Write header to file

for ih in xrange(shape(hgt)[0]):
    # For each phase of eruption height
    source = Sources(nLayers,name[ih])       # Initialize new source
    print(release_start)

    source.t0 = [release_start[ih] for i in xrange(nLayers)]
    source.t1 = [release_end[ih]   for i in xrange(nLayers)]
    type(hgt[ih])
    layers=[(hgt[ih]-summitHeight)*(i+1)/nLayers+summitHeight for i in ilayer]
    source.calcLayers(layers,mass[ih],p)

    source.fPrint(fout)             # Write source param to file

    
fout.close()

#dummy_source.fPrint(a)
#a.close()
#level=range(10)
#p = poisson(level,3)
#print "Poisson x:"
#print level
#print "Poisson probability density:"
#print p
