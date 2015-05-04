#!/usr/env python

import pprocess  # parallel computing module
import os  # os utilitites
import time  # keep track of time
import matplotlib.pyplot as plt
import numpy as np  # numeric python for array manipulation
from astropy.io import fits  # FITS manipulating library


def z2n(freqs, time, harm=1):
    '''
    z2n(freqs,time[, harm])

    Description
    -----------
    A python implementation of the Rayleigh Z_n^2 method to calculate the power
    spectrum of a time series in a given range of frequencies

    Caculates a periodogram, using Fourrier-analysis trough
    a method using the statistical variable Z^2_n wich has a
    probability density funcion equal to that of a X^2 (chi-squared)
    with 2n degrees of freedom

    Input
    ----------
    freqs: an array with frequencies in units of 1/s

    time: an array with the time series where to find a period

    harm: [optional] harmonic of the Fourrier analysis use higher harmonics to
    low signals.

    Returns
    ----------
      Z2n: An array with the powerspectrum of the time series

    '''
    N = len(time)
    Z2n = []
    for ni in freqs:
        aux = 0
        for k in xrange(harm):
            Phi = (ni*time) % 1
            arg = (k+1)*Phi*2.0*np.pi
            phicos = np.cos(arg)
            phisin = np.sin(arg)
            aux = aux + (phicos.sum())**2 + (phisin.sum())**2
        Z2n.append((2.0/N)*aux)

    return Z2n


def makefits(outptname='output.fits', *cols):
    '''
    makefits(outptname, col1[, col2, col3, ...])

    Description
    ----
    Creates a fits file with the provided columns information

    Input
    ----
    outptname : a string to use as name of the output file

    cols : one or more lists in the format:
        col1 : [array, 'name-of-the-column', 'format', 'unit']

    Returns
    ----
    Boolean True

    Output
    ----
    creates the file <outptname> on the current directory

    '''
    columns = []
    for i in xrange(len(cols)):
        columns.append(fits.Column(name=cols[i][1], format=cols[i][2],
                                   unit=cols[i][3], array=cols[i][0]))
    tbhdu = fits.TableHDU.from_columns(columns)
    tbhdu.writeto(outptname)
    print '\n Created file {0} \n'.format(outptname)

    return True


#input file informations and variables atributions---------
inptname = str(raw_input('Input file (with extension): '))
outptname = os.path.splitext(inptname)[0]+'_z2n_output.fits'
inpt = fits.open(inptname)
times = inpt[1].data.field('TIME')
inpt.close()

interval = float(times.max()-times.min())
startf = 1.0/interval

print "The start frequency is: ", startf
query = str(raw_input("Change the start frequency? (y/n): "))
if (query == 'y') or (query == 'Y'):
    startf = float(raw_input('Enter the start frequency: '))
else:
    pass

endf = float(raw_input('Enter the last frequency: '))

over = float(raw_input('Enter oversample factor: '))
fact = 1.0/over
deltaf = fact/interval

print "The frequency step will be: ", deltaf
query2 = str(raw_input('Change frequency step? (y/n): '))
if (query2 == 'y') or (query2 == 'Y'):
    deltaf = float(raw_input('Enter the frequency interval: '))
else:
    pass

freqs = np.arange(startf, endf, deltaf)

#harm = int(raw_input('The Harmonic to be considered:'))
harm = 1
#----------------- The parallelism starts here ------------------------
nproc = int(raw_input('Enter the number of processor to use: '))
if nproc < 1:
    nproc = 1  # default number of cpus = 1

freqlist = np.array_split(freqs, nproc)

results = pprocess.Map(limit=nproc, reuse=1)
parallel_z2n = results.manage(pprocess.MakeReusable(z2n))

print "\n Calculating with ", nproc, " processor(s)\n"

tic = time.time()
[parallel_z2n(somefreqs, times, harm) for somefreqs in freqlist]

z2n = []
for result in results:
    for value in result:
        z2n.append(value)

print 'time = {0}'.format(time.time() - tic)

plt.plot(freqs, z2n)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power')
plt.title(inptname)
plt.savefig('z2n_plot.pdf'.format(inpt), bbox_inches='tight', format='pdf', orientation='landscape')
plt.show()

#create and write the output.fits file
col1 = [freqs, 'frequency', 'E', 'Hz']
col2 = [z2n, 'z2nPower', 'E', 'Power']
makefits(outptname, col1, col2)
