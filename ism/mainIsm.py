
# MAIN FUNCTION TO CALL THE ISM MODULE

from ism.src.ism import ism

# Directory - this is the common directory for the execution of the E2E, all modules
auxdir = r'C:\Users\Kolarov\Desktop\EarthObservation\Codes\S3\test_eodp\auxiliary'
indir = r"C:\\Users\Kolarov\Desktop\EarthObservation\Codes\S3\EODP-TS-ISM\input\gradient_alt100_act150"
# indir = r"C:\\Users\Kolarov\Desktop\EarthObservation\Codes\S3\EODP-TS-E2E\sgm_out" # for E2E
outdir = r"C:\\Users\Kolarov\Desktop\EarthObservation\Codes\S3\EODP-TS-ISM\myoutput"
# outdir = r"C:\\Users\Kolarov\Desktop\EarthObservation\Codes\S3\EODP-TS-ISM\myoutput_E2E" # for E2E

# Initialise the ISM
myIsm = ism(auxdir, indir, outdir)
myIsm.processModule()

