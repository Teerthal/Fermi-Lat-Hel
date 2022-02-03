# Generate MC samples

# Usage: Read comments in bcut.py, qcalc.py for file structure
# This script takes two arguments: period and galactic latitude cut
# Example: python mc.py jan14z100 80
# generates the files in mc80.tar.gz and stores them under
# ./montecarlo/jan14z100


import numpy as np
import astropy.io.fits as pyfits #Replacing pyfits with astropy
import sys, os

# number of monte carlo samples
nummc = 10000

# deg to rad
degrad = np.pi / 180.

# obtain the files with the event data. We only need them to have the counts

#datadir = 'sourcecuts/' + sys.argv[1] + '/'
datadir = '/media/teerthal/Repo 2/LAT/'+sys.argv[1]+'/data_outputs/sourcecuts/'
file_list = []

# galactic latitude cut
galbcut = np.int(sys.argv[2])

for file in [doc for doc in os.listdir(datadir) if doc.startswith('nosource')]:
    file_list.append(file)

file_list = np.sort(file_list)


# print file_list

# return a matrix of unit vectors pointing towards every event
def evdirections(edata):
    ecb = np.cos(edata.field('b') * degrad)
    esb = np.sin(edata.field('b') * degrad)
    ecl = np.cos(edata.field('l') * degrad)
    esl = np.sin(edata.field('l') * degrad)
    return np.column_stack((ecl * ecb, esl * ecb, esb))


def mcsamplenorth(nevents, bmin, sourcemask):
    zmin = np.sin(bmin * degrad)
    i = 0
    mcdirections = np.zeros((nevents, 3))
    while i < nevents:
        # phi is uniformly distributed over 0, 2*pi
        phi = np.random.uniform(0., 2. * np.pi)
        # cos (theta) is uniformly distributed
        # we generate uniform numbers between zmin=sin(bmin) and 1
        # add random sign to sample both north and southern hemispheres
        # z=np.random.uniform(zmin,1)*np.random.choice([-1,1])
        z = np.random.uniform(zmin, 1)
        sqz = np.sqrt(1 - z ** 2)
        mcdirections[i] = [sqz * np.cos(phi), sqz * np.sin(phi), z]
        # we now have a uniformly distributed point in the sphere
        # we keep it if it is not too close to a source
        if np.dot(sourcemask, mcdirections[i]).max() < msize:
            i += 1  # keep direction if the largest value of cos(source,direction)
            # which is the closest to any source is below mcsize=3deg
    return mcdirections


def mcsamplesouth(nevents, bmin, sourcemask):
    zmin = np.sin(bmin * degrad)
    i = 0
    mcdirections = np.zeros((nevents, 3))
    while i < nevents:
        # phi is uniformly distributed over 0, 2*pi
        phi = np.random.uniform(0., 2. * np.pi)
        # cos (theta) is uniformly distributed
        # we generate uniform numbers between zmin=sin(bmin) and 1
        # add random sign to sample both north and southern hemispheres
        z = np.random.uniform(zmin, 1)
        sqz = np.sqrt(1 - z ** 2)
        mcdirections[i] = [sqz * np.cos(phi), sqz * np.sin(phi), -z]
        # we now have a uniformly distributed point in the sphere
        # we keep it if it is not too close to a source
        if np.dot(sourcemask, mcdirections[i]).max() < msize:
            i += 1  # keep direction if the largest value of cos(source,direction)
            # which is the closest to any source is below mcsize=3deg
    return mcdirections


# This is the angular region that we are going to consider 1-20 deg.
angstep = 1
minangle = 1
maxangle = 20

region = np.linspace(minangle, maxangle, int((maxangle - minangle + 1.) / angstep))

# the number of possible e1, e2 choices is len(file_list)-1
nume12 = len(file_list) - 1

# Source catalogue. We will avoid 3deg region around each source

# 1st HE source catalogue
#fermicatalogue = '/srv/data/fermi/sources/1fhel/gll_psch_v07.fit'
fermicatalogue='/media/teerthal/Repo 2/LAT/sources/gll_psch_v07.fit'

# 2nd source catalogue
# fermicatalogue='/srv/data/fermi/sources/2yrps/gll_psc_v08.fit'

fsource = pyfits.open(fermicatalogue)

# Size of mask around each source. This is the cosine.
msize = np.cos(1.5 * degrad)

# Create an array of unit vectors pointing in the direction of the sources
# can't use evdirections because the column is named glat, glon not l,b

# We only need sources above 48deg, since the closest gamma-ray is
# above 50deg. We separate north and south sources as well

srcdata = fsource[1].data
index = np.where(srcdata.field('glat') > 45.)
cb = np.cos(srcdata.field('glat')[index] * degrad)
sb = np.sin(srcdata.field('glat')[index] * degrad)
cl = np.cos(srcdata.field('glon')[index] * degrad)
sl = np.sin(srcdata.field('glon')[index] * degrad)

sourcenorth = np.column_stack((cl * cb, sl * cb, sb))

index = np.where(srcdata.field('glat') < -45.)
cb = np.cos(srcdata.field('glat')[index] * degrad)
sb = np.sin(srcdata.field('glat')[index] * degrad)
cl = np.cos(srcdata.field('glon')[index] * degrad)
sl = np.sin(srcdata.field('glon')[index] * degrad)

sourcesouth = np.column_stack((cl * cb, sl * cb, sb))

fsource.close()

fe3 = pyfits.open(datadir + file_list[len(file_list) - 1])
# number of E3 events above our cut of b
ne3n = np.shape(np.where(fe3[1].data.field('b') > galbcut))[1]
ne3s = np.shape(np.where(fe3[1].data.field('b') < -galbcut))[1]
fe3.close()

# we loop over the different combinations of E1 and E2
# for each one we generate mc montecarlos
# we save the values of q and delta_q

for i in range(nume12):
    for j in range(i + 1, nume12):

        # we will save the q, qnorth and qsouth values for each montecarlo here
        qmc = np.zeros((nummc, 3, len(region)))

        # obtain number of E1,E2 events above galbcut + maxangle(Region)
        fe1 = pyfits.open(datadir + file_list[i])
        ne1n = np.shape(np.where(fe1[1].data.field('b') > galbcut - maxangle))[1]
        ne1s = np.shape(np.where(fe1[1].data.field('b') < -galbcut + maxangle))[1]
        fe2 = pyfits.open(datadir + file_list[j])
        ne2n = np.shape(np.where(fe2[1].data.field('b') > galbcut - maxangle))[1]
        ne2s = np.shape(np.where(fe2[1].data.field('b') < -galbcut + maxangle))[1]
        fe1.close()
        fe2.close()
        # we have now the number of E1, E2 and E3 events. We do the MC
        for k in range(nummc):
            e3vecn = mcsamplenorth(ne3n, galbcut, sourcenorth)
            e3vecs = mcsamplesouth(ne3s, galbcut, sourcesouth)
            e1vecn = mcsamplenorth(ne1n, galbcut - maxangle, sourcenorth)
            e1vecs = mcsamplesouth(ne1s, galbcut - maxangle, sourcesouth)
            e2vecn = mcsamplenorth(ne2n, galbcut - maxangle, sourcenorth)
            e2vecs = mcsamplesouth(ne2s, galbcut - maxangle, sourcesouth)
            # matrix with all the scalar products of e3 and e1 or e2
            dote1e3n = np.dot(e3vecn, e1vecn.transpose())
            dote1e3s = np.dot(e3vecs, e1vecs.transpose())
            dote2e3n = np.dot(e3vecn, e2vecn.transpose())
            dote2e3s = np.dot(e3vecs, e2vecs.transpose())
            qval = np.zeros(len(region))
            qnorth = np.zeros(len(region))
            qsouth = np.zeros(len(region))
            for l in range(len(region)):  # start from larger region is faster
                rsize = np.cos(region[len(region) - 1 - l] * degrad)
                # mask pairs that are too far away
                e1inregion_n = dote1e3n > rsize
                e1inregion_s = dote1e3s > rsize
                e2inregion_n = dote2e3n > rsize
                e2inregion_s = dote2e3s > rsize

                eta1n = np.zeros((e3vecn.shape[0], 3))  # N_e3 values of eta1
                eta1s = np.zeros((e3vecs.shape[0], 3))  # N_e3 values of eta1
                eta2n = np.zeros((e3vecn.shape[0], 3))  # N_e3 values of eta1
                eta2s = np.zeros((e3vecs.shape[0], 3))  # N_e3 values of eta1

                eta1 = np.zeros((e3vecn.shape[0] + e3vecs.shape[0], 3))
                eta2 = np.zeros((e3vecn.shape[0] + e3vecs.shape[0], 3))

                for t in range(len(eta1n)):
                    if (e1vecn[e1inregion_n[t]].shape[0] == 0 or
                            e2vecn[e2inregion_n[t]].shape[0] == 0):
                        eta1n[t] = np.zeros(3)
                        eta2n[t] = np.zeros(3)
                    else:
                        eta1n[t] = np.average(e1vecn[e1inregion_n[t]], axis=0)
                        eta2n[t] = np.average(e2vecn[e2inregion_n[t]], axis=0)

                # calculate cross product
                eta1c2n = np.cross(eta1n, eta2n)  # array of N_e3 vectors

                # dot with eta3 and sum, then divide by n3. store in qval array
                eta1c2d3n = np.diag(np.dot(e3vecn, eta1c2n.transpose()))
                qnorth[len(region) - 1 - l] = eta1c2d3n.mean()

                for t in range(len(eta1s)):
                    if (e1vecs[e1inregion_s[t]].shape[0] == 0 or
                            e2vecs[e2inregion_s[t]].shape[0] == 0):
                        eta1s[t] = np.zeros(3)
                        eta2s[t] = np.zeros(3)
                    else:
                        eta1s[t] = np.average(e1vecs[e1inregion_s[t]], axis=0)
                        eta2s[t] = np.average(e2vecs[e2inregion_s[t]], axis=0)

                # calculate cross product
                eta1c2s = np.cross(eta1s, eta2s)  # array of N_e3 vectors

                # dot with eta3 and sum, then divide by n3. store in qval array
                eta1c2d3s = np.diag(np.dot(e3vecs, eta1c2s.transpose()))
                qsouth[len(region) - 1 - l] = eta1c2d3s.mean()

                eta1 = np.concatenate((eta1n, eta1s), axis=0)
                eta2 = np.concatenate((eta2n, eta2s), axis=0)

                eta1c2 = np.cross(eta1, eta2)
                e3vec = np.concatenate((e3vecn, e3vecs), axis=0)
                eta1c2d3 = np.diag(np.dot(e3vec, eta1c2.transpose()))
                qval[len(region) - 1 - l] = eta1c2d3.mean()
            # Multiply everything by 10^6
            qval *= 1.e6
            qnorth *= 1.e6
            qsouth *= 1.e6

            # We store the results in qmc matrix that we'll write to disc
            qmc[k] = [qval, qnorth, qsouth]

        np.save('/media/teerthal/Repo 2/LAT/'+sys.argv[1]+'/data_outputs/montecarlo/' #northsouthbgt' + str(galbcut)
                + str((i + 1) * 10) + '_' + str((j + 1) * 10) + 'GeV_nummc'
                + str(nummc), qmc)

