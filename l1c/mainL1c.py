
# MAIN FUNCTION TO CALL THE L1C MODULE

from l1c.src.l1c import l1c

# Directory - this is the common directory for the execution of the E2E, all modules
auxdir = r'C:\\Users\\Aless\\Desktop\\EARTH\\repo\\auxiliary'
# GM dir + L1B dir
indir = r'C:\\Users\\SHARED\\EODP-TS-E2E\\gm_out,C:\\Users\\SHARED\\EODP-TS-E2E\\output_2'
outdir = r'C:\\Users\\SHARED\\EODP-TS-E2E\\output_3'

# Initialise the ISM
myL1c = l1c(auxdir, indir, outdir)
myL1c.processModule()
