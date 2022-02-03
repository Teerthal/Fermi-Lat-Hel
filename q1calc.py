# Calculate Q' statistic from files with sources masked out

# Create plots for data above 70 and 80deg gal lat

# this script takes two arguments, period and galactic latitude cut
# For example: python qcalc.py jan14z100 80
# produces the plot jan14z100bgt80.pdf and the files qarrays.tar.gz

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
index = np.where(np.abs(e3data.field('b')) > galbcut)
e3vec70 = e3vec[index[0], :]

# This is the angular region that we are going to consider 1-20 deg.
angstep = 1
minangle = 1
maxangle = 20

region = np.linspace(minangle, maxangle, int((maxangle - minangle + 1.) / angstep))

# the number of possible e1, e2 files is len(file_list)-1
nume12 = len(file_list) - 1

import matplotlib
import matplotlib.pyplot as mp

def qcalc():
    matplotlib.use('pdf')
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
            dote1e3 = np.dot(e3vec70, e1vec.transpose())
            dote2e3 = np.dot(e3vec70, e2vec.transpose())
            qval = np.zeros(len(region))
            qerr = np.zeros(len(region))
            for k in range(len(region)):  # start from the larger region is faster
                msize = np.cos(region[len(region) - 1 - k] * degrad)
                # mask pairs that are too far away
                e1inregion = dote1e3 > msize
                e2inregion = dote2e3 > msize

                eta1 = np.zeros((e3vec70.shape[0], 3))  # N_e3 values of eta1
                eta2 = np.zeros((e3vec70.shape[0], 3))  # N_e3 values of eta1

                for t in range(len(eta1)):
                    if (e1vec[e1inregion[t]].shape[0] == 0 or
                            e2vec[e2inregion[t]].shape[0] == 0):
                        eta1[t] = np.zeros(3)
                        eta2[t] = np.zeros(3)
                    else:
                        eta1[t] = np.average(e1vec[e1inregion[t]], axis=0)
                        eta2[t] = np.average(e2vec[e2inregion[t]], axis=0)

                # calculate cross product
                eta1c2 = np.cross(eta1, eta2)  # this is an array of N_e3 vectors

                # dot with eta3 and sum, then divide by n3 and store in qval array
                eta1c2d3 = np.diag(np.dot(e3vec70, eta1c2.transpose()))
                qval[len(region) - 1 - k] = eta1c2d3.mean()
                qerr[len(region) - 1 - k] = eta1c2d3.std()
            # Divide the std in qerr by sqrt(N_e3)
            qerr /= np.sqrt(e3vec70.shape[0])
            # Multiply everything by 10^6
            qval *= 1.e6
            qerr *= 1.e6

            # We save the results to file in txt form
            qplot = np.column_stack((region, qval, qerr))
            np.savetxt('/media/teerthal/Repo 2/LAT/' + sys.argv[1] + '/data_outputs/qstat/' + sys.argv[1] + 'bgt' + str(
                galbcut) + str((i + 1) * 10)
                       + '_' + str((j + 1) * 10) + 'GeV.txt', qplot)

            # We now have all the values of q for this particular choice of E1,E2
            # Lets plot them
            ax[min(i, 1), spcol % 3].set_xlim((0, 21))
            ax[min(i, 1), spcol % 3].plot([0, 21], [0, 0], 'k-')
            ax[min(i, 1), spcol % 3].errorbar(region, qval, yerr=qerr, fmt='o', ecolor='g', )
            ax[min(i, 1), spcol % 3].set_title(r'$(E_1,E_2) = $' + '(' + str((i + 1) * 10) +
                                               ', ' + str((j + 1) * 10) + ')')

            spcol += 1

    # Fine-tune figure
    mp.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
    fig.suptitle(sys.argv[1] + ' Gal latidude > ' + str(galbcut), fontsize=22)
    fig.savefig('/media/teerthal/Repo 2/LAT/' + sys.argv[1] + '/data_outputs/qstat/' + sys.argv[1] + 'bgt' + str(
        galbcut) + '.pdf')

    return []

#qcalc()

def qcal(e1vec, e2vec, e3vec):
    # matrix with all the scalar products of e3 and e1 or e2
    dote1e3 = np.dot(e3vec, e1vec.transpose())
    dote2e3 = np.dot(e3vec, e2vec.transpose())
    qval = np.zeros(len(region))
    qerr = np.zeros(len(region))
    for k in range(len(region)):  # start from the larger region is faster
        msize = np.cos(region[len(region) - 1 - k] * degrad)
        # mask pairs that are too far away
        e1inregion = dote1e3 > msize
        e2inregion = dote2e3 > msize

        eta1 = np.zeros((e3vec70.shape[0], 3))  # N_e3 values of eta1
        eta2 = np.zeros((e3vec70.shape[0], 3))  # N_e3 values of eta1

        for t in range(len(eta1)):
            if (e1vec[e1inregion[t]].shape[0] == 0 or
                    e2vec[e2inregion[t]].shape[0] == 0):
                eta1[t] = np.zeros(3)
                eta2[t] = np.zeros(3)
            else:
                eta1[t] = np.average(e1vec[e1inregion[t]], axis=0)
                eta2[t] = np.average(e2vec[e2inregion[t]], axis=0)

        # calculate cross product
        eta1c2 = np.cross(eta1, eta2)  # this is an array of N_e3 vectors

        # dot with eta3 and sum, then divide by n3 and store in qval array
        eta1c2d3 = np.diag(np.dot(e3vec, eta1c2.transpose()))
        qval[len(region) - 1 - k] = eta1c2d3.mean()
        qerr[len(region) - 1 - k] = eta1c2d3.std()
    # Divide the std in qerr by sqrt(N_e3)
    qerr /= np.sqrt(e3vec.shape[0])
    # Multiply everything by 10^6
    qval *= 1.e6
    qerr *= 1.e6

    return [qval, qerr]

def newq(e1vec, e2vec, e3vec):
    # matrix with all the scalar products of e3 and e1 or e2
    dote1e3 = np.dot(e3vec, e1vec.transpose())
    dote2e3 = np.dot(e3vec, e2vec.transpose())
    qval = np.zeros(len(region))
    qerr = np.zeros(len(region))
    for l in range(len(region)):  # start from larger region is faster
        rsize = np.cos(np.deg2rad(region[len(region) - 1 - l]))
        # mask pairs that are too far away
        e1inregion = dote1e3 > rsize
        e2inregion = dote2e3 > rsize

        # now we need to take care of the theta function for each E3

        q_ser = np.zeros(e3vec.shape[0])   #Storing series for computing error in q
        for i in range(e3vec.shape[0]):  # loop over each E3
            e1in = e1vec[e1inregion[i]]  # possible E1 within region around E3
            e2in = e2vec[e2inregion[i]]

            if (e1in.shape[0] > 0 and e2in.shape[0] > 0):  # zero if no E1 or E2
                # build m1 and m2. z direction is given by e3 vector
                m1 = e1in - np.dot(np.dot(e1in, e3vec[i]).reshape((e1in.shape[0], 1)), e3vec[i].reshape((1, 3)))
                m2 = e2in - np.dot(np.dot(e2in, e3vec[i]).reshape((e2in.shape[0], 1)), e3vec[i].reshape((1, 3)))

                # the theta term is a NE1 x NE2 matrix
                # we define it as a boolean matrix: true when the scalar
                # product is positive
                # each row gives us the locations of the E2 directions
                # that are relevant, so we can sum those as if it where eta_2

                theta = np.dot(m1, m2.transpose()) > 0

                eta2 = np.zeros((e1in.shape[0], 3))

                for j in range(len(eta2)):
                    eta2[j] = e2in[theta[j]].sum(axis=0)
                    # if(len(e2in[theta[j]])>0):
                    #   eta2[j] = np.average(e2in[theta[j]],axis=0)

                # now we cross each E1 with the corresponding eta2, dot it with
                # E3 and sum
                if (np.sum(theta) > 0):
                    term = np.sum(np.dot(np.cross(e1in, eta2), e3vec[i]), axis=0) / np.sum(theta)
                    qval[len(region) - 1 - l] += term
                else:
                    term = 0
                    qval[len(region) - 1 - l] += term
                # qval[l]/=(e1in.shape[0] * e2in.shape[0]) # divide by N1 N2

            q_ser[i] = term
        # once we finish the loop over all E3, we divide by NE3
        qval[len(region) - 1 - l] /= e3vec.shape[0]
        qerr[len(region) - 1 - l] = q_ser.std()/np.sqrt(e3vec.shape[0])
    qval *= 1.e6
    qerr *= 1.e6
    return [qval, qerr]

def rel_plot():

    ####Plot new and old q values on same tiles######

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

            #Calling Q and Q' values

            q_mat = qcal(e1vec, e2vec, e3vec70)
            q_val = q_mat[0]
            q_err = q_mat[1]

            q_mat_1 = newq(e1vec, e2vec, e3vec70)
            qval_1 = q_mat_1[0]
            qerr_1 = q_mat_1[1]

            # Lets plot them
            ax[min(i, 1), spcol % 3].set_xlim((0, 21))
            ax[min(i, 1), spcol % 3].plot([0, 21], [0, 0], 'k-')
            ax[min(i, 1), spcol % 3].errorbar(region, q_val, yerr=q_err, fmt='ro', ecolor='r', )
            #markers, caps, bars = \
            #    ax[min(i, 1), spcol % 3].errorbar(region, qval_1, yerr=qerr_1, fmt='go', ecolor='g', elinewidth = 2, alpha = 0.5)
            ax[min(i, 1), spcol % 3].set_title(r'$(E_1,E_2) = $' + '(' + str((i + 1) * 10) +
                                               ', ' + str((j + 1) * 10) + ')')
            #[bar.set_alpha(0.5) for bar in bars]
            #[cap.set_alpha(0.5) for cap in caps]

            spcol += 1

    # Fine-tune figure
    mp.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
    #fig.suptitle(sys.argv[1] + ' Gal latidude > ' + str(galbcut), fontsize=22)
    fig.savefig('/media/teerthal/Repo 2/LAT/' + sys.argv[1] + '/data_outputs/qstat/' + sys.argv[1] + 'bgt' + str(
        galbcut) + '.pdf')
    mp.show()
    return

rel_plot()