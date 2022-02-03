# Calculate Q statistic from files with sources masked out

# Create plots for data above 70 and 80deg gal lat

import numpy as np
import astropy.io.fits as pyfits #Replacing pyfits with astropy
import sys, os

# deg to rad
degrad = np.pi / 180.

# obtain the files with the event data

datadir = '/media/teerthal/Repo 2/LAT/' + sys.argv[1] + '/data_outputs/sourcecuts/'

file_list = []

# galactic latitude cut
galbcut = np.int(sys.argv[2])

for file in [doc for doc in os.listdir(datadir) if doc.startswith('nosource')]:
    file_list.append(file)

file_list = np.sort(file_list)

print
file_list


# return a matrix of unit vectors pointing towards every event
def evdirections(edata):
    ecb = np.cos(edata.field('b') * degrad)
    esb = np.sin(edata.field('b') * degrad)
    ecl = np.cos(edata.field('l') * degrad)
    esl = np.sin(edata.field('l') * degrad)
    return np.column_stack((ecl * ecb, esl * ecb, esb))


fe3 = pyfits.open(datadir + file_list[len(file_list) - 1])
e3data = fe3[1].data
e3vec = evdirections(e3data)
# we only want e3 events above 70 and 80 deg gal latitude
index = np.where(e3data.field('b') > galbcut)
e3vec70north = e3vec[index[0], :]
index = np.where(e3data.field('b') < -galbcut)
e3vec70south = e3vec[index[0], :]

# This is the angular region that we are going to consider 1-20 deg.
angstep = 1
minangle = 1
maxangle = 20

region = np.linspace(minangle, maxangle, int((maxangle - minangle + 1.) / angstep))

# the number of possible e1, e2 files is len(file_list)-1
nume12 = len(file_list) - 1

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as mp
fig, ax = mp.subplots(nrows=2, ncols=3, sharex='col')
spcol = 0

for i in range(nume12):
    for j in range(i + 1, nume12):
        fe1 = pyfits.open(datadir + file_list[i])
        fe2 = pyfits.open(datadir + file_list[j])
        e1data = fe1[1].data
        e1vec = evdirections(e1data)
        e2data = fe2[1].data
        e2vec = evdirections(e2data)
        # matrix with all the scalar products of e3 and e1 or e2
        dote1e3north = np.dot(e3vec70north, e1vec.transpose())
        dote2e3north = np.dot(e3vec70north, e2vec.transpose())
        qvalnorth = np.zeros(len(region))
        qerrnorth = np.zeros(len(region))
        dote1e3south = np.dot(e3vec70south, e1vec.transpose())
        dote2e3south = np.dot(e3vec70south, e2vec.transpose())
        qvalsouth = np.zeros(len(region))
        qerrsouth = np.zeros(len(region))
        for k in range(len(region)):  # start from the larger region is faster
            msize = np.cos(region[len(region) - 1 - k] * degrad)
            # mask pairs that are too far away
            e1inregionnorth = dote1e3north > msize
            e2inregionnorth = dote2e3north > msize

            eta1north = np.zeros((e3vec70north.shape[0], 3))  # N_e3 values of eta1
            eta2north = np.zeros((e3vec70north.shape[0], 3))  # N_e3 values of eta1

            for t in range(len(eta1north)):
                if (e1vec[e1inregionnorth[t]].shape[0] == 0 or
                        e2vec[e2inregionnorth[t]].shape[0] == 0):
                    eta1north[t] = np.zeros(3)
                    eta2north[t] = np.zeros(3)
                else:
                    eta1north[t] = np.average(e1vec[e1inregionnorth[t]], axis=0)
                    eta2north[t] = np.average(e2vec[e2inregionnorth[t]], axis=0)

            # calculate cross product
            eta1c2north = np.cross(eta1north, eta2north)  # this is an array of N_e3 vectors

            # dot with eta3 and sum, then divide by n3 and store in qval array
            eta1c2d3north = np.diag(np.dot(e3vec70north, eta1c2north.transpose()))
            qvalnorth[len(region) - 1 - k] = eta1c2d3north.mean()
            qerrnorth[len(region) - 1 - k] = eta1c2d3north.std()

            # repeat for the southern hemisphere
            e1inregionsouth = dote1e3south > msize
            e2inregionsouth = dote2e3south > msize

            eta1south = np.zeros((e3vec70south.shape[0], 3))  # N_e3 values of eta1
            eta2south = np.zeros((e3vec70south.shape[0], 3))  # N_e3 values of eta1

            for t in range(len(eta1south)):
                if (e1vec[e1inregionsouth[t]].shape[0] == 0 or
                        e2vec[e2inregionsouth[t]].shape[0] == 0):
                    eta1south[t] = np.zeros(3)
                    eta2south[t] = np.zeros(3)
                else:
                    eta1south[t] = np.average(e1vec[e1inregionsouth[t]], axis=0)
                    eta2south[t] = np.average(e2vec[e2inregionsouth[t]], axis=0)

            # calculate cross product
            eta1c2south = np.cross(eta1south, eta2south)  # this is an array of N_e3 vectors

            # dot with eta3 and sum, then divide by n3 and store in qval array
            eta1c2d3south = np.diag(np.dot(e3vec70south, eta1c2south.transpose()))
            qvalsouth[len(region) - 1 - k] = eta1c2d3south.mean()
            qerrsouth[len(region) - 1 - k] = eta1c2d3south.std()
        # Divide the std in qerr by sqrt(N_e3)
        qerrnorth /= np.sqrt(e3vec70north.shape[0])
        qerrsouth /= np.sqrt(e3vec70south.shape[0])
        # Multiply everything by 10^6
        qvalnorth *= 1.e6
        qvalsouth *= 1.e6
        qerrnorth *= 1.e6
        qerrsouth *= 1.e6

        # We save the results to file in txt form
        qplot = np.column_stack((region, qvalnorth, qerrnorth, qvalsouth, qerrsouth))
        np.savetxt('/media/teerthal/Repo 2/LAT/' + sys.argv[1] + '/data_outputs/qstat/northsouth' + sys.argv[1] + 'bgt' + str(galbcut) + str((i + 1) * 10)
                   + '_' + str((j + 1) * 10) + 'GeV.txt', qplot)
        np.save('/media/teerthal/Repo 2/LAT/' + sys.argv[1] + '/data_outputs/qstat/northsouth' + sys.argv[1] + 'bgt' + str(galbcut) + str((i + 1) * 10)
                   + '_' + str((j + 1) * 10) + 'GeV', qplot)

        # We now have all the values of q for this particular choice of E1,E2
        # Lets plot them
        ax[min(i, 1), spcol % 3].set_xlim((0, 21))
        ax[min(i, 1), spcol % 3].plot([0, 21], [0, 0], 'k-')
        ax[min(i, 1), spcol % 3].errorbar(region, qvalnorth, yerr=qerrnorth, fmt='go', ecolor='g', )
        ax[min(i, 1), spcol % 3].errorbar(region, qvalsouth, yerr=qerrsouth, fmt='ro', ecolor='r', )
        ax[min(i, 1), spcol % 3].set_title(r'$(E_1,E_2) = $' + '(' + str((i + 1) * 10) +
                                           ', ' + str((j + 1) * 10) + ')')

        spcol += 1

# Fine-tune figure
mp.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
fig.suptitle(sys.argv[1] + ' Gal latidude > ' + str(galbcut), fontsize=22)
fig.savefig('/media/teerthal/Repo 2/LAT/' + sys.argv[1] + '/data_outputs/qstat/northsouth' + sys.argv[1] + 'bgt' + str(galbcut) + '.pdf')

