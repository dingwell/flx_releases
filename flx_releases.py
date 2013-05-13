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
        

    def calcLayers(self,hgts,mass,scaleFactors):
        #TODO: check arguments
        # input hgt in [m]
        # sets up the source layers based on total column height
        # hgts is a list of vertical layer tops
        # p is an array with scale factors (from probability density function)
        # p will be normalized so that sum(p)=1
        # We also require the global variable summitHeight to determine the
        # bottom of the lowermost layer.
        
        # Normalize scaleFactors
        sfm = scaleFactors/sum(scaleFactors)        # for Mass
        sfr = [math.sqrt(i) for i in scaleFactors]  # for Length (box side)
        sfr = sfr/sum(sfr)
        # The reason we use different scale factors for mass and length:
        # We want concentration to be about constant throughout the eruption
        # column, concentration in each layer is a function of mass and radius:
        #   C = const*mass/radius^2     [ m/(pi*r^2*h) ]
        # (layer thickness, h, is constant)
        # Therefore, maximum concentration would be at the smallest r.
        # By using the sqrt of the sf for r we get something more realistic

        self.z1[:]  = [h*1000 for h in hgts[:]]       # Top of layers
        self.z0[1:] = [h*1000 for h in hgts[:-1]]     # Bottom of layers
        self.z0[0]  = summitHeight*1000  # Bottom of lowest layer

        for ih in xrange(shape(hgts)[0]):
            # For each layer (slice) of the source
            h = hgts[ih]*1000   # We want height in m
            self.mass[ih] = mass*sfm[ih]
            dx,dy = self.metres2degrees(h*sfr[ih])
            self.x0[ih] = centerLon-dx/2
            self.x1[ih] = centerLon+dx/2
            self.y0[ih] = centerLat-dy/2
            self.y1[ih] = centerLat+dy/2

            # Get duration of eruption phase
            #diff = [i-j for i,j in zip(self.t1,self.t0)]    # duration as timedelta
            #dt   = [d.days*24+d.seconds/3600 for d in diff] # duration in hours
            diff = self.t1[ih]-self.t0[ih]
            dt   = float(diff.days*24)+float(diff.seconds)/3600.0
            print dt,diff.seconds
            if dt <= 0.0:
                print "t0,t1:",self.t1[ih],"\t",self.t0[ih]
                print "dt:",dt
                raise Exception("Duration of phase is too short, make sure"\
                                "t0 and t1 are properly defined")

            self.npart[ih] = npart/self.n*dt # N particles to be released
            

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
