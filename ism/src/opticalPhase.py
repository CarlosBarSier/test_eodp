from ism.src.initIsm import initIsm
from math import pi
from ism.src.mtf import mtf
from numpy.fft import fftshift, ifftshift, ifft2, fft2
import numpy as np
from common.io.writeToa import writeToa
from common.io.readIsrf import readIsrf
from scipy.interpolate import interp1d, interp2d
from common.plot.plotMat2D import plotMat2D
from common.plot.plotF import plotF
from scipy.signal import convolve2d
from common.src.auxFunc import getIndexBand

class opticalPhase(initIsm):

    def __init__(self, auxdir, indir, outdir):
        super().__init__(auxdir, indir, outdir)

    def compute(self, sgm_toa, sgm_wv, band):
        """
        The optical phase is in charge of simulating the radiance
        to irradiance conversion, the spatial filter (PSF)
        and the spectral filter (ISRF).
        :return: TOA image in irradiances [mW/m2/nm],
                    with spatial and spectral filter
        """
        self.logger.info("EODP-ALG-ISM-1000: Optical stage")

        # Calculation and application of the ISRF
        # -------------------------------------------------------------------------------
        self.logger.info("EODP-ALG-ISM-1010: Spectral modelling. ISRF")
        toa = self.spectralIntegration(sgm_toa, sgm_wv, band)

        self.logger.debug("TOA [0,0] " +str(toa[0,0]) + " [e-]")

        if self.ismConfig.save_after_isrf:
            saveas_str = self.globalConfig.ism_toa_isrf + band
            writeToa(self.outdir, saveas_str, toa)

        # Radiance to Irradiance conversion
        # -------------------------------------------------------------------------------
        self.logger.info("EODP-ALG-ISM-1020: Radiances to Irradiances")
        toa = self.rad2Irrad(toa,
                             self.ismConfig.D,
                             self.ismConfig.f,
                             self.ismConfig.Tr)

        self.logger.debug("TOA [0,0] " +str(toa[0,0]) + " [e-]")

        # Spatial filter
        # -------------------------------------------------------------------------------
        # Calculation and application of the system MTF
        self.logger.info("EODP-ALG-ISM-1030: Spatial modelling. PSF/MTF")
        myMtf = mtf(self.logger, self.outdir)
        Hsys = myMtf.system_mtf(toa.shape[0], toa.shape[1],
                                self.ismConfig.D, self.ismConfig.wv[getIndexBand(band)], self.ismConfig.f, self.ismConfig.pix_size,
                                self.ismConfig.kLF, self.ismConfig.wLF, self.ismConfig.kHF, self.ismConfig.wHF,
                                self.ismConfig.defocus, self.ismConfig.ksmear, self.ismConfig.kmotion,
                                self.outdir, band)

        # Apply system MTF
        toa = self.applySysMtf(toa, Hsys) # always calculated
        self.logger.debug("TOA [0,0] " +str(toa[0,0]) + " [e-]")



        # Write output TOA & plots
        # -------------------------------------------------------------------------------
        if self.ismConfig.save_optical_stage:
            saveas_str = self.globalConfig.ism_toa_optical + band

            writeToa(self.outdir, saveas_str, toa)

            title_str = 'TOA after the optical phase [mW/sr/m2]'
            xlabel_str='ACT'
            ylabel_str='ALT'
            plotMat2D(toa, title_str, xlabel_str, ylabel_str, self.outdir, saveas_str)

            idalt = int(toa.shape[0]/2)
            saveas_str = saveas_str + '_alt' + str(idalt)
            plotF([], toa[idalt,:], title_str, xlabel_str, ylabel_str, self.outdir, saveas_str)

            # --- DEBUG OPTICAL STAGE ---
            print("\n[DEBUG OpticalStage]")
            print("TOA shape:", toa.shape)
            print("TOA mean/std antes del MTF (estimado):", np.mean(sgm_toa), np.std(sgm_toa))
            print("TOA después del MTF:", np.mean(toa), np.std(toa))
            print("Hsys min/max:", np.min(Hsys), np.max(Hsys))
            print("Ruta de salida:", self.outdir)
            print("----------------------------\n")

        return toa

    def rad2Irrad(self, toa, D, f, Tr):
        """
        Radiance to Irradiance conversion
        :param toa: Input TOA image in radiances [mW/sr/m2]
        :param D: Pupil diameter [m]
        :param f: Focal length [m]
        :param Tr: Optical transmittance [-]
        :return: TOA image in irradiances [mW/m2]
        """

        toa = Tr*toa*(np.pi/4)*(D/f)**2
        return toa


    def applySysMtf(self, toa, Hsys):
        """
        Application of the system MTF to the TOA
        :param toa: Input TOA image in irradiances [mW/m2]
        :param Hsys: System MTF
        :return: TOA image in irradiances [mW/m2]
        """
        """
        Apply the system MTF to the TOA image (same dimensions).
        The MTF provided by mtf.system_mtf() is already centred.
        """
        i0, j0 = toa.shape[0] // 2, toa.shape[1] // 2
        dc = Hsys[i0, j0]
        if dc != 0:
            Hsys = Hsys / dc
        # Aplica SIN shifts (Hsys ya está centrado como lo generas)
        GE = fft2(toa)
        toa_ft = np.real(ifft2(GE * Hsys))
        return toa_ft

    def spectralIntegration(self, sgm_toa, sgm_wv, band):
        """
        Integration with the ISRF to retrieve one band
        :param sgm_toa: Spectrally oversampled TOA cube 3D in irradiances [mW/m2]
        :param sgm_wv: wavelengths of the input TOA cube
        :param band: band
        :return: TOA image 2D in radiances [mW/m2]
        """
        # wv in [um]
        isrf, wv_isrf_um = readIsrf(self.auxdir + self.ismConfig.isrffile, band)
        wv_isrf_nm = wv_isrf_um * 1000.0

        # Normaliza el ISRF por su INTEGRAL (no por la suma)
        # -> área = 1 en el eje de wv de trabajo (nm aquí)
        area = np.trapz(isrf, wv_isrf_nm)
        isrf_n = isrf / area if area != 0 else isrf

        toa = np.zeros((sgm_toa.shape[0], sgm_toa.shape[1]))
        for ialt in range(sgm_toa.shape[0]):
            for iact in range(sgm_toa.shape[1]):
                # interpola la columna espectral del SGM a la rejilla del ISRF
                cs = interp1d(sgm_wv, sgm_toa[ialt, iact, :],
                              fill_value=(0.0, 0.0), bounds_error=False)
                toa_i = cs(wv_isrf_nm)  # TOA(λ) en la rejilla del ISRF
                # integra TOA(λ)*ISRF_n(λ) dλ
                toa[ialt, iact] = np.trapz(toa_i * isrf_n, wv_isrf_nm)

        return toa
