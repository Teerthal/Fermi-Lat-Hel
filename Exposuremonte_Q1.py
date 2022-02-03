# Calculate Q statistic from files with sources masked out

# Create plots for data above 70 and 80deg gal lat

import numpy as np
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
import sys, os

# number of monte carlo samples
nummc = 10

# deg to rad
degrad = np.pi / 180.

# obtain the files with the event data. We only need them to have the counts

datadir = '/home/tpatel28/Fermi/' + sys.argv[1] + '/sourcecuts/'
#datadir='/media/teerthal/Repo 2/LAT/'+sys.argv[1]+'/data_outputs/sourcecuts/'


file_list = []

# galactic latitude cut
galbcut = np.int(sys.argv[2])

for file in [doc for doc in os.listdir(datadir) if doc.startswith('nosource')]:
    file_list.append(file)

file_list = np.sort(file_list)

# print file_list

# obtain the exposure files.
expdir = '/home/tpatel28/Fermi/' + sys.argv[1] + '/exp_maps/'
#expdir='/media/teerthal/Repo 2/LAT/'+sys.argv[1]+'/data_outputs/exp_maps/'
file_exp = []

for file in [doc for doc in os.listdir(expdir) if doc.startswith('all_sky_')]:
    file_exp.append(file)

file_exp = np.sort(file_exp)

# we build the wcs structure by hand

wcs = pywcs.WCS(naxis=2)
wcs.wcs.crpix = [900.5, 450.5]
wcs.wcs.cdelt = np.array([-0.2, 0.2])
wcs.wcs.crval = [0, 0]
wcs.wcs.ctype = ["GLON-AIT", "GLAT-AIT"]


# return a matrix of unit vectors pointing towards every event
def evdirections(edata):
    ecb = np.cos(edata.field('b') * degrad)
    esb = np.sin(edata.field('b') * degrad)
    ecl = np.cos(edata.field('l') * degrad)
    esl = np.sin(edata.field('l') * degrad)
    return np.column_stack((ecl * ecb, esl * ecb, esb))


# now the mc samples take an exposure array too
def mcsample(nevents, bmin, sourcemask, exposure):
    zmin = np.sin(bmin * degrad)
    i = 0
    mcdirections = np.zeros((nevents, 3))
    while i < nevents:
        # phi is uniformly distributed over 0, 2*pi
        phi = np.random.uniform(0., 2. * np.pi)
        # cos (theta) is uniformly distributed
        # we generate uniform numbers between zmin=sin(bmin) and 1
        # add random sign to sample both north and southern hemispheres
        z = np.random.uniform(zmin, 1) * np.random.choice([-1, 1])
        sqz = np.sqrt(1 - z ** 2)
        mcdirections[i] = [sqz * np.cos(phi), sqz * np.sin(phi), z]
        # calculate angular coordinates in degrees for pixel mapping
        phi /= degrad  # this is the galactic longitude
        bgal = np.arctan(z / sqz) / degrad  # galactic latitude in degrees

        # pixcoord=wcs.wcs_sky2pix(np.array([[phi,bgal]]),1)
        # ### wcs_sky2pix is phased out for wcs_world2pix

        pixcoord = wcs.wcs_world2pix(np.array([[phi, bgal]]), 1)

        # note that the pixel coordinates are swapped in the numpy array

        # we now have a uniformly distributed point in the sphere
        # we keep it if it is not too close to a source and
        # we modulate it by the exposure which is normalized between 0 and 1
        if (np.dot(sourcemask, mcdirections[i]).max() < msize and
                exposure[int(pixcoord[0, 1]), int(pixcoord[0, 0])] >
                np.random.uniform(0, 1)):
            i += 1  # keep direction if the largest value of cos(source,direction)
            # which is the closest to any source is below mcsize=3deg

            ####Hotfixed by Teerthal above since the pixcoordinates seemed not to have integer values######
            ####Need to check if world2pix is working as intended###########

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
fermicatalogue = '/home/tpatel28/Fermi/sources/gll_psch_v07.fit'
#fermicatalogue='/media/teerthal/Repo 2/LAT/sources/gll_psch_v07.fit'
# 2nd source catalogue
# fermicatalogue='/home/tpatel28/Fermi/sources/gll_psc_v08.fit'

# 3rd HE source catalogue
fermicatalogue='/home/tpatel28/Fermi/sources/gll_psch_v23.fit'

fsource = pyfits.open(fermicatalogue)

# Size of mask around each source. This is the cosine.
msize = np.cos(1.5 * degrad)

# Create an array of unit vectors pointing in the direction of the sources
# can't use evdirections because the column is named glat, glon not l,b
srcdata = fsource[1].data
cb = np.cos(srcdata.field('glat') * degrad)
sb = np.sin(srcdata.field('glat') * degrad)
cl = np.cos(srcdata.field('glon') * degrad)
sl = np.sin(srcdata.field('glon') * degrad)

sourcevec = np.column_stack((cl * cb, sl * cb, sb))

fsource.close()

fe3 = pyfits.open(datadir + file_list[len(file_list) - 1])
# number of E3 events above our cut of b
ne3 = np.shape(np.where(np.abs(fe3[1].data.field('b')) > galbcut))[1]
fe3.close()

fexp3 = pyfits.open(expdir + file_exp[len(file_exp) - 1])
dataexp3 = fexp3[0].data.mean(axis=0)  # we average over energy bins
dataexp3 /= dataexp3.max()  # normalize to 1
fexp3.close()

# we loop over the different combinations of E1 and E2
# for each one we generate mc montecarlos
# we save the values of q and delta_q
print(nume12)
for i in range(nume12):

    for j in range(i + 1, nume12):

        # we will save the q and delta_q values for each montecarlo here
        qmc = np.zeros((nummc, 2, len(region)))

        # obtain number of E1,E2 events above galbcut + maxangle(Region)
        fe1 = pyfits.open(datadir + file_list[i])
        ne1 = np.shape(np.where(np.abs(fe1[1].data.field('b')) >
                                galbcut - maxangle))[1]
        fe2 = pyfits.open(datadir + file_list[j])
        ne2 = np.shape(np.where(np.abs(fe2[1].data.field('b')) >
                                galbcut - maxangle))[1]
        fe1.close()
        fe2.close()

        # get the exposures
        fexp1 = pyfits.open(expdir + file_exp[i])
        dataexp1 = fexp1[0].data.mean(axis=0)  # we average over energy bins
        dataexp1 /= dataexp1.max()  # normalize to 1
        fexp1.close()
        fexp2 = pyfits.open(expdir + file_exp[j])
        dataexp2 = fexp2[0].data.mean(axis=0)  # we average over energy bins
        dataexp2 /= dataexp2.max()  # normalize to 1
        fexp2.close()

        # we have now the number of E1, E2 and E3 events. We do the MC
        for k in range(nummc):
            e3vec = mcsample(ne3, galbcut, sourcevec, dataexp3)
            e1vec = mcsample(ne1, galbcut - maxangle, sourcevec, dataexp1)
            e2vec = mcsample(ne2, galbcut - maxangle, sourcevec, dataexp2)
            # matrix with all the scalar products of e3 and e1 or e2
            dote1e3 = np.dot(e3vec, e1vec.transpose())
            dote2e3 = np.dot(e3vec, e2vec.transpose())
            qval = np.zeros(len(region))
            qerr = np.zeros(len(region))
            for l in range(len(region)):  # start from larger region is faster
                rsize = np.cos(region[len(region) - 1 - l] * degrad)
                # mask pairs that are too far away
                e1inregion = dote1e3 > rsize
                e2inregion = dote2e3 > rsize

                q_ser = np.zeros(e3vec.shape[0])  # Storing series for computing error in q
                for idx in range(e3vec.shape[0]):  # loop over each E3
                    e1in = e1vec[e1inregion[idx]]  # possible E1 within region around E3
                    e2in = e2vec[e2inregion[idx]]

                    if (e1in.shape[0] > 0 and e2in.shape[0] > 0):  # zero if no E1 or E2
                        # build m1 and m2. z direction is given by e3 vector
                        m1 = e1in - np.dot(np.dot(e1in, e3vec[idx]).reshape((e1in.shape[0], 1)), e3vec[idx].reshape((1, 3)))
                        m2 = e2in - np.dot(np.dot(e2in, e3vec[idx]).reshape((e2in.shape[0], 1)), e3vec[idx].reshape((1, 3)))

                        # the theta term is a NE1 x NE2 matrix
                        # we define it as a boolean matrix: true when the scalar
                        # product is positive
                        # each row gives us the locations of the E2 directions
                        # that are relevant, so we can sum those as if it where eta_2

                        theta = np.dot(m1, m2.transpose()) > 0

                        eta2 = np.zeros((e1in.shape[0], 3))

                        for jdx in range(len(eta2)):
                            eta2[jdx] = e2in[theta[jdx]].sum(axis=0)
                            # if(len(e2in[theta[j]])>0):
                            #   eta2[j] = np.average(e2in[theta[j]],axis=0)

                        # now we cross each E1 with the corresponding eta2, dot it with
                        # E3 and sum
                        if (np.sum(theta) > 0):
                            term = np.sum(np.dot(np.cross(e1in, eta2), e3vec[idx]), axis=0) / np.sum(theta)
                            qval[len(region) - 1 - l] += term
                        else:
                            term = 0
                            qval[len(region) - 1 - l] += term
                        # qval[l]/=(e1in.shape[0] * e2in.shape[0]) # divide by N1 N2

                    q_ser[idx] = term
                # once we finish the loop over all E3, we divide by NE3
                qval[len(region) - 1 - l] /= e3vec.shape[0]
                qerr[len(region) - 1 - l] = q_ser.std() / np.sqrt(e3vec.shape[0])

            qval *= 1.e6
            qerr *= 1.e6

            # We store the results in qmc matrix that we'll write to disc
            qmc[k] = [qval, qerr]

        np.save('/home/tpatel28/Fermi/' + sys.argv[1] + '/exp_mc_1_ser/exposurebgt' + str(galbcut)
                + str((i + 1) * 10) + '_' + str((j + 1) * 10) + 'GeV_nummc' + str(nummc), qmc)

        #np.save('/media/teerthal/Repo 2/LAT/' + sys.argv[1] + '/data_outputs/montecarlo/Exp_MC_1/exposurebgt' + str(galbcut)
         #       + str((i + 1) * 10) + '_' + str((j + 1) * 10) + 'GeV_nummc' + str(nummc), qmc)
