from math import pi
from config.ismConfig import ismConfig
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import j1
from numpy.matlib import repmat
from common.io.readMat import writeMat
from common.plot.plotMat2D import plotMat2D
from scipy.interpolate import interp2d
from numpy.fft import fftshift, ifft2
import os

class mtf:
    """
    Class MTF. Collects the analytical modelling of the different contributions
    for the system MTF
    """
    def __init__(self, logger, outdir):
        self.ismConfig = ismConfig()
        self.logger = logger
        self.outdir = outdir

    def system_mtf(self, nlines, ncolumns, D, lambd, focal, pix_size,
                   kLF, wLF, kHF, wHF, defocus, ksmear, kmotion, directory, band):
        """
        System MTF
        :param nlines: Lines of the TOA
        :param ncolumns: Columns of the TOA
        :param D: Telescope diameter [m]
        :param lambd: central wavelength of the band [m]
        :param focal: focal length [m]
        :param pix_size: pixel size in meters [m]
        :param kLF: Empirical coefficient for the aberrations MTF for low-frequency wavefront errors [-]
        :param wLF: RMS of low-frequency wavefront errors [m]
        :param kHF: Empirical coefficient for the aberrations MTF for high-frequency wavefront errors [-]
        :param wHF: RMS of high-frequency wavefront errors [m]
        :param defocus: Defocus coefficient (defocus/(f/N)). 0-2 low defocusing
        :param ksmear: Amplitude of low-frequency component for the motion smear MTF in ALT [pixels]
        :param kmotion: Amplitude of high-frequency component for the motion smear MTF in ALT and ACT
        :param directory: output directory
        :return: mtf
        """

        self.logger.info("Calculation of the System MTF")

        # Calculate the 2D relative frequencies
        self.logger.debug("Calculation of 2D relative frequencies")
        fn2D, fr2D, fnAct, fnAlt = self.freq2d(nlines, ncolumns, D, lambd, focal, pix_size)

        # Diffraction MTF
        self.logger.debug("Calculation of the diffraction MTF")
        Hdiff = self.mtfDiffract(fr2D)

        # Defocus
        Hdefoc = self.mtfDefocus(fr2D, defocus, focal, D)

        # WFE Aberrations
        Hwfe = self.mtfWfeAberrations(fr2D, lambd, kLF, wLF, kHF, wHF)

        # Detector
        Hdet  = self. mtfDetector(fn2D)

        # Smearing MTF
        Hsmear = self.mtfSmearing(fnAlt, ncolumns, ksmear)

        # Motion blur MTF
        Hmotion = self.mtfMotion(fn2D, kmotion)

        # Calculate the System MTF
        self.logger.debug("Calculation of the Sysmtem MTF by multiplying the different contributors")
        Hsys = Hdiff*Hdefoc*Hwfe*Hdet*Hsmear*Hmotion # REVISAR

        # Plot cuts ACT/ALT of the MTF
        self.plotMtf(Hdiff, Hdefoc, Hwfe, Hdet, Hsmear, Hmotion, Hsys, nlines, ncolumns, fnAct, fnAlt, directory, band)


        return Hsys

    def freq2d(self, nlines, ncolumns, D, lambd, focal, w):
        """
        Calculate the relative/normalised spatial frequencies (2D and 1D)
        """
        import numpy as np

        eps = 1e-6


        fstepAlt = 1.0 / (nlines * w)
        fstepAct = 1.0 / (ncolumns * w)

        # --- SIN np.arange PORQUE PROBLEMAS COMPATIBILIDAD
        idxAlt = np.array(list(range(nlines)), dtype=float)
        idxAct = np.array(list(range(ncolumns)), dtype=float)

        fAlt = (idxAlt - nlines / 2.0) * fstepAlt
        fAct = (idxAct - ncolumns / 2.0) * fstepAct

        fAlt[-1] = fAlt[-1] - eps
        fAct[-1] = fAct[-1] - eps

        fnAlt = fAlt * w
        fnAct = fAct * w

        fc = D / (lambd * focal)
        frAlt = fAlt / fc
        frAct = fAct / fc

        fnAltxx, fnActxx = np.meshgrid(fnAlt, fnAct, indexing='ij')
        frAltxx, frActxx = np.meshgrid(frAlt, frAct, indexing='ij')

        fn2D = np.sqrt(fnAltxx * fnAltxx + fnActxx * fnActxx)
        fr2D = np.sqrt(frAltxx * frAltxx + frActxx * frActxx)

        return fn2D, fr2D, fnAct, fnAlt

    #def freq2d(self,nlines, ncolumns, D, lambd, focal, w):
        """
        Calculate the relative frequencies 2D (for the diffraction MTF)
        :param nlines: Lines of the TOA
        :param ncolumns: Columns of the TOA
        :param D: Telescope diameter [m]
        :param lambd: central wavelength of the band [m]
        :param focal: focal length [m]
        :param w: pixel size in meters [m]
        :return fn2D: normalised frequencies 2D (f/(1/w))
        :return fr2D: relative frequencies 2D (f/(1/fc))
        :return fnAct: 1D normalised frequencies 2D ACT (f/(1/w))
        :return fnAlt: 1D normalised frequencies 2D ALT (f/(1/w))
        """

        fstepAlt = 1/nlines/w
        fstepAct = 1/ncolumns/w
        cutoff = D/lambd/focal
        eps = 1e-10
        fAlt = np.arrange(-1/(2*w),1/(2*w)-eps,fstepAlt)
        fAct = np.arrange(-1 / (2 * w), 1 / (2 * w) - eps, fstepAct)

        #1D vectors
        #frAlt normalized with cut-off
        #frAct
        #frAlt normalized with 1/w
        #frAct

        fnAct = fAct/(1/w)
        fnAlt = fAlt/(1/w)

        #2D frequency vectors
        [fnAltxx, fnActxx] = np.meshgrid(fnAlt, fnAct, indexing='ij')
        fn2D = np.sqrt(fnAltxx*fnAltxx + fnActxx*fnActxx)
        fr2D = fn2D*(1/w)/cutoff
        #2D fr

        return fn2D, fr2D, fnAct, fnAlt

    def mtfDiffract(self,fr2D):
        """
        Optics Diffraction MTF
        :param fr2D: 2D relative frequencies (f/fc), where fc is the optics cut-off frequency
        :return: diffraction MTF
        """
        Hdiff = np.zeros_like(fr2D)
        for i in range(fr2D.shape[0]):
            for j in range(fr2D.shape[1]):
                if fr2D[i, j] < 1:
                    Hdiff[i, j] = 2 / np.pi * (
                            np.arccos(fr2D[i,j]) - fr2D[i,j] * np.sqrt(1 - fr2D[i,j] ** 2)
                    )
                else:
                    Hdiff[i, j] = 0.0
        return Hdiff


    def mtfDefocus(self, fr2D, defocus, focal, D):
        """
        Defocus MTF
        :param fr2D: 2D relative frequencies (f/fc), where fc is the optics cut-off frequency
        :param defocus: Defocus coefficient (defocus/(f/N)). 0-2 low defocusing
        :param focal: focal length [m]
        :param D: Telescope diameter [m]
        :return: Defocus MTF
        """
        x = np.pi * defocus * fr2D * (1-fr2D)
        J1 = x/2 - x**3/16 + x**5/384 - x**7/18432
        Hdefoc = 2*J1/(x+1e-12)
        return Hdefoc

    def mtfWfeAberrations(self, fr2D, lambd, kLF, wLF, kHF, wHF):
        """
        Wavefront Error Aberrations MTF
        :param fr2D: 2D relative frequencies (f/fc), where fc is the optics cut-off frequency
        :param lambd: central wavelength of the band [m]
        :param kLF: Empirical coefficient for the aberrations MTF for low-frequency wavefront errors [-]
        :param wLF: RMS of low-frequency wavefront errors [m]
        :param kHF: Empirical coefficient for the aberrations MTF for high-frequency wavefront errors [-]
        :param wHF: RMS of high-frequency wavefront errors [m]
        :return: WFE Aberrations MTF
        """
        a = kLF*(wLF*wLF/lambd/lambd) + kHF*(wHF*wHF/lambd/lambd)
        Hwfe = np.exp(-fr2D * (1 - fr2D) * a)
        return Hwfe

    def mtfDetector(self,fn2D):
        """
        Detector MTF
        :param fnD: 2D normalised frequencies (f/(1/w))), where w is the pixel width
        :return: detector MTF
        """
        Hdet = np.abs(np.sin(np.pi*fn2D)/(np.pi*fn2D + 1e-12))
        return Hdet

    def mtfSmearing(self, fnAlt, ncolumns, ksmear):
        """
        Smearing MTF
        :param ncolumns: Size of the image ACT
        :param fnAlt: 1D normalised frequencies 2D ALT (f/(1/w))
        :param ksmear: Amplitude of low-frequency component for the motion smear MTF in ALT [pixels]
        :return: Smearing MTF
        """
        Hsmear = np.zeros((fnAlt.shape[0], ncolumns))
        row = np.sin(fnAlt * ksmear * np.pi) / (fnAlt * ksmear * np.pi + 1e-12)
        for j in range(ncolumns):
            Hsmear[:, j] = row
        return Hsmear

    def mtfMotion(self, fn2D, kmotion):
        """
        Motion blur MTF
        :param fnD: 2D normalised frequencies (f/(1/w))), where w is the pixel width
        :param kmotion: Amplitude of high-frequency component for the motion smear MTF in ALT and ACT
        :return: detector MTF
        """
        Hmotion = np.sin(fn2D*kmotion*np.pi)/(fn2D*kmotion*np.pi + 1e-12)
        return Hmotion

    def plotMtf(self,Hdiff, Hdefoc, Hwfe, Hdet, Hsmear, Hmotion, Hsys, nlines, ncolumns, fnAct, fnAlt, directory, band):
        """
        Plotting the system MTF and all of its contributors
        :param Hdiff: Diffraction MTF
        :param Hdefoc: Defocusing MTF
        :param Hwfe: Wavefront electronics MTF
        :param Hdet: Detector MTF
        :param Hsmear: Smearing MTF
        :param Hmotion: Motion blur MTF
        :param Hsys: System MTF
        :param nlines: Number of lines in the TOA
        :param ncolumns: Number of columns in the TOA
        :param fnAct: normalised frequencies in the ACT direction (f/(1/w))
        :param fnAlt: normalised frequencies in the ALT direction (f/(1/w))
        :param directory: output directory
        :param band: band
        :return: N/A
        """
        #TODO
        #REVISA AL FINAL!!!!!


