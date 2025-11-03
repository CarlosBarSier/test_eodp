
# MAIN FUNCTION TO CALL THE L1B MODULE

from l1b.src.l1b import l1b

# Directory - this is the common directory for the execution of the E2E, all modules
auxdir = r'C:\\Users\Kolarov\Desktop\EarthObservation\Codes\S3\test_eodp\auxiliary'
indir = r"C:\\Users\Kolarov\Desktop\EarthObservation\Codes\S3\EODP-TS-L1B\input"
# outdir = r"C:\\Users\Kolarov\Desktop\EarthObservation\Codes\S3\EODP-TS-L1B\myoutputs_eq" # for equalization
outdir = r"C:\\Users\Kolarov\Desktop\EarthObservation\Codes\S3\EODP-TS-L1B\myoutputs_noeq" # for no equalization

# Initialise the ISM
myL1b = l1b(auxdir, indir, outdir)
myL1b.processModule()
