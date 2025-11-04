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

        writeMat(self.outdir, "fn2D", fn2D)
        writeMat(self.outdir, "fr2D", fr2D)

        return fn2D, fr2D, fnAct, fnAlt

    def mtfDiffract(self,fr2D):
        """
        Optics Diffraction MTF
        :param fr2D: 2D relative frequencies (f/fc), where fc is the optics cut-off frequency
        :return: diffraction MTF
        """
        #Hdiff = np.zeros_like(fr2D)
        #for i in range(fr2D.shape[0]):
        #    for j in range(fr2D.shape[1]):
        #        if fr2D[i, j] < 1:
        #            Hdiff[i, j] = 2 / np.pi * (
        #                    np.arccos(fr2D[i,j]) - fr2D[i,j] * np.sqrt(1 - fr2D[i,j] ** 2)
        #            )
        #        else:
        #            Hdiff[i, j] = 0.0

        fr = np.clip(fr2D, 0, 1)
        Hdiff = (2 / np.pi) * (np.arccos(fr) - fr * np.sqrt(1 - fr ** 2))
        Hdiff[fr2D >= 1] = 0.0

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
        #x = np.pi*defocus*fr2D*(1-fr2D)
        #J1 = x/2 - x**3/16 + x**5/384 - x**7/18432
        #Hdefoc = 2*J1/(x+1e-12)
        #return Hdefoc

        x = np.pi * defocus * fr2D * (1 - fr2D)
        Hdefoc = np.ones_like(x)
        mask = x != 0
        Hdefoc[mask] = 2 * j1(x[mask]) / x[mask]
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
        Hwfe = np.exp(-fr2D*(1-fr2D)*a)
        return Hwfe

    def mtfDetector(self,fn2D):
        """
        Detector MTF
        :param fnD: 2D normalised frequencies (f/(1/w))), where w is the pixel width
        :return: detector MTF
        """
        x = np.pi * fn2D
        Hdet = np.ones_like(fn2D)
        non_zero = ~np.isclose(x, 0.0)
        Hdet[non_zero] = np.abs(np.sin(x[non_zero]) / x[non_zero])
        return Hdet

    def mtfSmearing(self, fnAlt, ncolumns, ksmear):
        """
        Smearing MTF
        :param ncolumns: Size of the image ACT
        :param fnAlt: 1D normalised frequencies 2D ALT (f/(1/w))
        :param ksmear: Amplitude of low-frequency component for the motion smear MTF in ALT [pixels]
        :return: Smearing MTF
        """
        smear_1d = np.sinc(ksmear * fnAlt)  # valor 1 en fnAlt==0
        Hsmear = np.transpose(np.tile(smear_1d, (ncolumns, 1)))
        return Hsmear

    def mtfMotion(self, fn2D, kmotion):
        """
        Motion blur MTF
        :param fnD: 2D normalised frequencies (f/(1/w))), where w is the pixel width
        :param kmotion: Amplitude of high-frequency component for the motion smear MTF in ALT and ACT
        :return: detector MTF
        """
        Hmotion = np.sinc(kmotion * fn2D)           # 1 en el origen
        return Hmotion

    #def plotMtf(self,Hdiff, Hdefoc, Hwfe, Hdet, Hsmear, Hmotion, Hsys, nlines, ncolumns, fnAct, fnAlt, directory, band):
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

    def plotMtf(self, Hdiff, Hdefoc, Hwfe, Hdet, Hsmear, Hmotion, Hsys,
                    nlines, ncolumns, fnAct, fnAlt, directory, band):
            """
            Guarda:
              1) Hsys 2D (mapa)      -> ism_mtf2d_<band>.png
              2) Cortes ACT y ALT    -> ism_mtf_cut_act_<band>.png / ism_mtf_cut_alt_<band>.png
                 (incluye cada contribución y el sistema)
              3) MTFs en .mat        -> ism_H<component>_<band>.mat y ism_Hsys_<band>.mat
            Los ejes de los cortes usan frecuencia normalizada (f * w), coherente con fnAct/fnAlt.
            """
            import os
            import numpy as np
            import matplotlib.pyplot as plt

            # --- 0) Definir nombres de salida
            base2d = f"ism_mtf2d_{band}"
            base_act = f"ism_mtf_cut_act_{band}"
            base_alt = f"ism_mtf_cut_alt_{band}"

            # --- 1) Guardar los arrays en .mat (por si el test/validador los lee)
            try:
                writeMat(os.path.join(directory, f"ism_Hdiff_{band}.mat"), "Hdiff", Hdiff)
                writeMat(os.path.join(directory, f"ism_Hdefoc_{band}.mat"), "Hdefoc", Hdefoc)
                writeMat(os.path.join(directory, f"ism_Hwfe_{band}.mat"), "Hwfe", Hwfe)
                writeMat(os.path.join(directory, f"ism_Hdet_{band}.mat"), "Hdet", Hdet)
                writeMat(os.path.join(directory, f"ism_Hsmear_{band}.mat"), "Hsmear", Hsmear)
                writeMat(os.path.join(directory, f"ism_Hmotion_{band}.mat"), "Hmotion", Hmotion)
                writeMat(os.path.join(directory, f"ism_Hsys_{band}.mat"), "Hsys", Hsys)
            except Exception as e:
                self.logger.warning(f"No se pudieron guardar los .mat de MTF: {e}")

            # --- 2) Mapa 2D del MTF del sistema
            try:
                title_str = "System MTF"
                xlabel_str = "ACT (normalised frequency)"
                ylabel_str = "ALT (normalised frequency)"
                # Usamos plotMat2D utilitario del proyecto
                plotMat2D(Hsys, title_str, xlabel_str, ylabel_str, directory, base2d)
            except Exception as e:
                self.logger.warning(f"No se pudo generar el mapa 2D del MTF: {e}")

            # --- 3) Cortes 1D por el centro
            iact_c = int(ncolumns // 2)  # corte ALT: columna central
            ialt_c = int(nlines // 2)  # corte ACT: fila central

            # Cortes por componentes (para superponer curvas)
            components = {
                "Diffraction": Hdiff,
                "Defocus": Hdefoc,
                "WFE": Hwfe,
                "Detector": Hdet,
                "Smearing": Hsmear,
                "Motion": Hmotion,
                "System": Hsys
            }

            # --- 3a) Corte ACT (fila central): usar solo frecuencias >= 0
            iact_c = int(ncolumns // 2)
            ialt_c = int(nlines // 2)

            # índices de la mitad positiva (incluye 0)
            idx_act_pos = fnAct >= 0
            x_act = fnAct[idx_act_pos]

            plt.figure(figsize=(8, 5))
            for name, H in components.items():
                y_full = H[ialt_c, :]
                y_act = y_full[idx_act_pos]
                plt.plot(x_act, y_act, label=name)
            plt.title("System MTF, ACT cut (central ALT)")
            plt.xlabel("Spatial frequency (normalised, f·w)")
            plt.ylabel("MTF")
            plt.grid(True, alpha=0.3)
            plt.legend(loc="best", fontsize=8)
            plt.tight_layout()
            plt.savefig(os.path.join(directory, f"ism_mtf_cut_act_{band}.png"), dpi=150)
            plt.close()

            # --- 3b) Corte ALT (columna central): usar solo frecuencias >= 0
            idx_alt_pos = fnAlt >= 0
            x_alt = fnAlt[idx_alt_pos]

            plt.figure(figsize=(8, 5))
            for name, H in components.items():
                y_full = H[:, iact_c]
                y_alt = y_full[idx_alt_pos]
                plt.plot(x_alt, y_alt, label=name)
            plt.title("System MTF, ALT cut (central ACT)")
            plt.xlabel("Spatial frequency (normalised, f·w)")
            plt.ylabel("MTF")
            plt.grid(True, alpha=0.3)
            plt.legend(loc="best", fontsize=8)
            plt.tight_layout()
            plt.savefig(os.path.join(directory, f"ism_mtf_cut_alt_{band}.png"), dpi=150)
            plt.close()

            # --- Guardar valores MTF en Nyquist ---
            try:
                # Cálculo de Nyquist (frecuencia normalizada 0.5)
                target_fn = 0.5

                # Buscar el índice de Nyquist más cercano en los ejes ACT y ALT
                idx_act_ny = int(np.argmin(np.abs(fnAct - target_fn)))
                idx_alt_ny = int(np.argmin(np.abs(fnAlt - target_fn)))

                # Valores correspondientes
                mtf_act_val = float(Hsys[ialt_c, idx_act_ny])  # dirección ACT
                mtf_alt_val = float(Hsys[idx_alt_ny, iact_c])  # dirección ALT

                # Archivo de salida
                output_path = self.outdir + '/mtf_nyquist.txt'

                # Escribir encabezado solo una vez
                if band == 'VNIR-0':
                    with open(output_path, 'w') as file:
                        file.write('=== MTF at Nyquist (f·w = 0.5) ===\n\n')

                # Añadir línea con los valores de esta banda
                with open(output_path, 'a') as file:
                    file.write(f'{band}\n')
                    file.write(f'Act_MTF = {mtf_act_val:.6f}\n')
                    file.write(f'Alt_MTF = {mtf_alt_val:.6f}\n\n')

            except Exception as e:
                self.logger.warning(f'No se pudo escribir mtf_nyquist.txt: {e}')

