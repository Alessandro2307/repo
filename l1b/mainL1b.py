
# MAIN FUNCTION TO CALL THE L1B MODULE

from l1b.src.l1b import l1b

# Directory - this is the common directory for the execution of the E2E, all modules
auxdir = r'C:\\Users\\Aless\\Desktop\\EARTH\\repo\\auxiliary'
indir = r'C:\\Users\\SHARED\\EODP-TS-E2E\\output'
outdir = r"C:\\Users\\SHARED\\EODP-TS-E2E\\output_2"

# Initialise the ISM
myL1b = l1b(auxdir, indir, outdir)
myL1b.processModule()
