# Calculate Q statistic from files with sources masked out 

# We generate the MCs for all the panels at once. E.g. all 20GeV photons
# are the same for all pannels.


import numpy as np
#import astropy.io.fits as pyfits
import astropy.wcs as pywcs 
#import sys,os


# This is the angular region that we are going to consider 1-30 deg.
angstep=1
minangle=1
maxangle=30

region=np.linspace(minangle,maxangle,int((maxangle-minangle+1.)/angstep))

# random number generation initialization

prng = np.random.RandomState()

twopi=2*np.pi  #in appears as a bound on mc rng routines

# return a matrix of unit vectors pointing towards every event
def evdirections(edata):
    ecb = np.cos(np.deg2rad(edata.field('b')))
    esb = np.sin(np.deg2rad(edata.field('b')))
    ecl = np.cos(np.deg2rad(edata.field('l')))
    esl = np.sin(np.deg2rad(edata.field('l')))
    return np.column_stack((ecl*ecb,esl*ecb,esb))

# function that calculates q from the event directions
def calcq(e1vec,e2vec,e3vec):
    # matrix with all the scalar products of e3 and e1 or e2
    dote1e3=np.dot(e3vec,e1vec.transpose())
    dote2e3=np.dot(e3vec,e2vec.transpose())
    qval=np.zeros(len(region))
    qerr=np.zeros(len(region))
    for l in range(len(region)): #start from larger region is faster
        rsize = np.cos(np.deg2rad(region[len(region)-1-l]))
        # mask pairs that are too far away
        e1inregion = dote1e3 > rsize
        e2inregion = dote2e3 > rsize
            
        eta1=np.zeros((e3vec.shape[0],3)) # N_e3 values of eta1
        eta2=np.zeros((e3vec.shape[0],3)) # N_e3 values of eta1
            
        for t in range(len(eta1)):
            if (e1vec[e1inregion[t]].shape[0]==0 or 
                    e2vec[e2inregion[t]].shape[0]==0):
                eta1[t]=np.zeros(3)
                eta2[t]=np.zeros(3)
            else:
                eta1[t]=np.average(e1vec[e1inregion[t]],axis=0)
                eta2[t]=np.average(e2vec[e2inregion[t]],axis=0)

        #calculate cross product
        eta1c2=np.cross(eta1, eta2) # this is an array of N_e3 vectors
            
        #dot with eta3 and sum, then divide by n3. store in qval array
        eta1c2d3=np.diag(np.dot(e3vec,eta1c2.transpose()))
        qval[len(region)-1-l]=eta1c2d3.mean()
        qerr[len(region)-1-l]=eta1c2d3.std()
            	
# Divide the std in qerr by sqrt(N_e3)
    qerr/=np.sqrt(e3vec.shape[0])
    # Multiply everything by 10^6
    qval*=1.e6
    qerr*=1.e6

    return [qval,qerr]

#### New Q with the theta function


def newq(e1vec,e2vec,e3vec):
    # matrix with all the scalar products of e3 and e1 or e2
    dote1e3=np.dot(e3vec,e1vec.transpose())
    dote2e3=np.dot(e3vec,e2vec.transpose())
    qval=np.zeros(len(region))
    qerr=np.zeros(len(region))
    for l in range(len(region)): #start from larger region is faster
        rsize = np.cos(np.deg2rad(region[len(region)-1-l]))
        # mask pairs that are too far away
        e1inregion = dote1e3 > rsize
        e2inregion = dote2e3 > rsize
            
        # now we need to take care of the theta function for each E3
        
        for i in range(e3vec.shape[0]): # loop over each E3
            e1in = e1vec[e1inregion[i]] # possible E1 within region around E3
            e2in = e2vec[e2inregion[i]]

            if (e1in.shape[0]>0 and e2in.shape[0]>0): # zero if no E1 or E2
                # build m1 and m2. z direction is given by e3 vector
                m1 = e1in - np.dot( np.dot(e1in,e3vec[i]).reshape((e1in.shape[0],1)), e3vec[i].reshape((1,3)))
                m2 = e2in - np.dot( np.dot(e2in,e3vec[i]).reshape((e2in.shape[0],1)), e3vec[i].reshape((1,3)))

                # the theta term is a NE1 x NE2 matrix 
                # we define it as a boolean matrix: true when the scalar
                # product is positive
                # each row gives us the locations of the E2 directions
                # that are relevant, so we can sum those as if it where eta_2

                theta = np.dot(m1,m2.transpose()) > 0  

                eta2 = np.zeros((e1in.shape[0],3))

                for j in range(len(eta2)):
                    eta2[j] = e2in[theta[j]].sum(axis=0)
                    #if(len(e2in[theta[j]])>0):
                     #   eta2[j] = np.average(e2in[theta[j]],axis=0)
                        

                # now we cross each E1 with the corresponding eta2, dot it with
                # E3 and sum
                if(np.sum(theta)>0):
                    qval[len(region)-1-l]+=np.sum(np.dot( np.cross(e1in,eta2), e3vec[i]),axis=0)/np.sum(theta)
                else:
                    qval[len(region)-1-l]+=0
                #qval[l]/=(e1in.shape[0] * e2in.shape[0]) # divide by N1 N2

        # once we finish the loop over all E3, we divide by NE3
        qval[len(region)-1-l]/=e3vec.shape[0]
        
    qval*=1.e6
    qerr*=1.e6
    return [qval,qerr]

def calcqFD(e1vec,e2vec,e3vec):

        if(len(e1vec)*len(e2vec)*len(e3vec)==0):
            return [np.zeros(len(region)),np.zeros(len(region))]
        # matrix with all the scalar products of e3 and e1 or e2
        dote1e3=np.dot(e3vec,e1vec.transpose())
        dote2e3=np.dot(e3vec,e2vec.transpose())
        qval=np.zeros(len(region))
        qerr=np.zeros(len(region))
        for l in range(len(region)): #start from larger gv.region is faster
            rsize = np.cos(np.deg2rad(region[len(region)-1-l]))
            # mask pairs that are too far away
            e1inregion = dote1e3 > rsize
            e2inregion = dote2e3 > rsize
            
            eta1=np.zeros((e3vec.shape[0],3)) # N_e3 values of eta1
            eta2=np.zeros((e3vec.shape[0],3)) # N_e3 values of eta2
            Se2ce1=np.zeros((e3vec.shape[0],3)) # N_e3 values of sum e2 cross eta1
	    
            for t in range(len(eta1)):
                e1=e1vec[e1inregion[t]]
                e2=e2vec[e2inregion[t]]
                if (e1.shape[0]==0 or
                        e2.shape[0]==0):
                    eta1[t]=np.zeros(3)
                    eta2[t]=np.zeros(3)
                    Se2ce1[t]=np.zeros(3)
                else:
                    de2e3=dote2e3[t][e2inregion[t]]
                    de1e3=dote1e3[t][e1inregion[t]] 
                    m1=e1-np.outer(de1e3,e3vec[t])
                    m2=e2-np.outer(de2e3,e3vec[t])
                    inreg=np.dot(m2,m1.transpose())>0
                    e1s=np.zeros((len(inreg),3))
                    for intofIR in range(len(inreg)):
                        if(len(e1[inreg[intofIR]])>0):
                            e1s[intofIR]=np.sum(e1[inreg[intofIR]],axis=0)
                    Se2ce1[t]=np.zeros(3)
                    if(np.sum(inreg)>0):
                        Se2ce1[t]=-np.sum(np.cross(e2,e1s),axis=0)/np.sum(inreg) #the negative is to have the same sign.

            eta1c2d3=np.diag(np.dot(e3vec,Se2ce1.transpose()))
            qval[len(region)-1-l]=eta1c2d3.mean()
            qerr[len(region)-1-l]=eta1c2d3.std()
            	
# Divide the std in qerr by sqrt(N_e3)
        qerr/=np.sqrt(e3vec.shape[0])
        # Multiply everything by 10^6
        qval*=1.e6
        qerr*=1.e6
        return [qval,qerr]


####
####
## 3. Read the exposure files and create wcs 
####
####


def buildwcs(header): #reads the information from the fits file, returns wcs
    wcstemp=pywcs.WCS(header)
    wcs = pywcs.WCS(naxis=2)
    wcs.wcs.crpix=wcstemp.wcs.crpix[:2]
    wcs.wcs.cdelt=wcstemp.wcs.cdelt[:2]
    wcs.wcs.crval=wcstemp.wcs.crval[:2]
    wcs.wcs.ctype=[wcstemp.wcs.ctype[0],wcstemp.wcs.ctype[1]]
    return wcs



####
####
## 4. Generate the MCs. We do N&S separately
####
####

# function to generate montecarlos. Takes an exposure file and a 
# possible -1 for the south. It also takes a wcs structure

def mcsample(nevents,bmin,sourcemask,exposure,wcs,ns):

    # ns can be 1 or -1 depending on N or South

    zmin=np.sin(np.deg2rad(bmin))
    i=0
    mcdirections=np.zeros((nevents,3))
    while i < nevents:
        # phi is uniformly distributed over 0, 2*pi
        phi=prng.uniform(0.,twopi)
        # cos (theta) is uniformly distributed
        # we generate uniform numbers between zmin=sin(bmin) and 1
        # add random sign to sample both north and southern hemispheres
        z=prng.uniform(zmin,1)*ns # if ns=-1 we are in the South
        sqz=np.sqrt(1-z**2)
        mcdirections[i]=[sqz * np.cos(phi), sqz*np.sin(phi),z]
        # we now have a uniformly distributed point in the sphere
        # we keep it if it is not too close to a source and
        # we modulate it by the exposure which is normalized between 0 and 1
        # the exposure modulation into a sepparate if, so that we
        # only call the coordinate-pixel transf. if not excluded by sources
        if (np.dot(sourcemask,mcdirections[i]).max() < msize):
            # candidate direction is not blocked by a source,
            # the largest value of cos(source,direction)
            # which is the closest to any source is below mcsize=3deg,
            # proceed to modulate by exposure

            # calculate angular coordinates in degrees for pixel mapping
            phi=np.rad2deg(phi) # this is the galactic longitude
#           if(phi > 180.):  # it does not seem to be required in the end
#	    	phi-=360. # put phi in the range (-180,180) for wcs
            bgal=np.rad2deg(np.arctan2(z,sqz)) #galactic latitude in degrees
            pixcoord=wcs.wcs_world2pix(np.array([[phi,bgal]]),0)
            # note that the pixel coordinates are swapped in the numpy array
            if (exposure[int(pixcoord[0,1]),int(pixcoord[0,0])] > 
                    prng.uniform(0,1)):
                i+=1 #keep direction
    return mcdirections


def mcsamplenoexp(nevents,bmin,sourcemask,ns):

    # ns can be 1 or -1 depending on N or South

    zmin=np.sin(np.deg2rad(bmin))
    
    i=0
    mcdirections=np.zeros((nevents,3))
    while i < nevents:
        # phi is uniformly distributed over 0, 2*pi
        phi=prng.uniform(0.,twopi)
        # cos (theta) is uniformly distributed
        # we generate uniform numbers between zmin=sin(bmin) and 1
        # add random sign to sample both north and southern hemispheres
        z=prng.uniform(zmin,1)*ns # if ns=-1 we are in the South
        sqz=np.sqrt(1-z**2)
        mcdirections[i]=[sqz * np.cos(phi), sqz*np.sin(phi),z]
       # we now have a uniformly distributed point in the sphere
        # we keep it if it is not too close to a source and
        if (np.dot(sourcemask,mcdirections[i]).max() < msize):
            i+=1 #keep direction if the largest value of cos(source,direction)
            # which is the closest to any source is below mcsize=3deg
    return mcdirections



def mcsampleuniform(nevents,bmin,ns):

    # ns can be 1 or -1 depending on N or South

    zmin=np.sin(np.deg2rad(bmin))
    
    i=0
    mcdirections=np.zeros((nevents,3))
    while i < nevents:
        # phi is uniformly distributed over 0, 2*pi
        phi=prng.uniform(0.,twopi)
        # cos (theta) is uniformly distributed
        # we generate uniform numbers between zmin=sin(bmin) and 1
        # add random sign to sample both north and southern hemispheres
        z=prng.uniform(zmin,1)*ns # if ns=-1 we are in the South
        sqz=np.sqrt(1-z**2)
        mcdirections[i]=[sqz * np.cos(phi), sqz*np.sin(phi),z]
    return mcdirections



