import os
import sys
import re
import blibrary
import argparse
import warnings
import linecache
import multiprocessing
import itertools

os.system("taskset -p 0xfffff %d" % os.getpid())
os.environ['MKL_NUM_THREADS'] = '1' 
os.environ['NUMEXPR_NUM_THREADS'] = '1' 
os.environ['OMP_NUM_THREADS'] = '1' 

sys.stdout.flush()

parser = argparse.ArgumentParser(description='This will binodals.')
parser.add_argument('--addr-stable',   metavar='astable',   dest='astable', type=str, action='store', help='enter address of stable islands.')
parser.add_argument('--addr-unstable', metavar='aunstable', dest='austble', type=str, action='store', help='enter address of unstable islands.')
parser.add_argument('--addr-crits',    metavar='acrits',    dest='acrits',  type=str, action='store', help='enter address of critical points.')
parser.add_argument('--addr-binodals', metavar='abins',     dest='abins',   type=str, action='store', help='enter address of binodals points.')
parser.add_argument('--mesh',          metavar='mesh',      dest='mesh',    type=str, action='store', help="Name of mesh pkl.")
parser.add_argument('--nproc',         metavar='nproc',     dest='nproc',   type=int, action='store', help="Call these many processes.")
parser.add_argument('--binodal-pkl',   metavar='bpkl',      dest='bpkl',    type=str, action='store', help="Name of binodal pkl.", default=None)
parser.add_argument('--no-rtw',         dest='nrtw',       action='store_true',  help="Dont print out the runtime warning.", default=False)
args = parser.parse_args()

#########################################
def custom_warning_format(message, category, filename, lineno, line=None):
    line = linecache.getline(filename, lineno).strip()
    if args.nrtw:
        return f""
    else:
        return f"There is a RunTimeWarning taking place on line {lineno}.\n"

warnings.formatwarning = custom_warning_format
#########################################

if __name__=="__main__":

    l_stable   = os.listdir(args.astable)
    l_unstable = os.listdir(args.austble)
    l_crits    = os.listdir(args.acrits)
    l_binodals = os.listdir(args.abins)
    nproc      = args.nproc

    files = []

    completed_calculations = []
    for file in l_binodals:
        print(f"file = {file}")
        chips = (re.search("chips_(-?\d+(?:\.\d+)?)", file)).group(1)
        chipc = (re.search("chipc_(-?\d+(?:\.\d+)?)", file)).group(1)
        chisc = (re.search("chisc_(-?\d+(?:\.\d+)?)", file)).group(1)
        vs    = (re.search("vs_(-?\d+(?:\.\d+)?)", file)).group(1)
        vc    = (re.search("vc_(-?\d+(?:\.\d+)?)", file)).group(1)
        vp    = (re.search("vp_(-?\d+(?:\.\d+)?)", file)).group(1)
        completed_calculations.append([float(chips), float(chipc), float(chisc), float(vs), float(vc), float(vp)])

    same_count = 0
    for file in l_unstable:
        to_del = re.search("TODEL", file)
        if to_del:
            # print(f"Ignoring {file}...")
            continue
        else:
            flist = [file]
            chips = (re.search("chips_(-?\d+(?:\.\d+)?)", file)).group(1)
            chipc = (re.search("chipc_(-?\d+(?:\.\d+)?)", file)).group(1)
            chisc = (re.search("chisc_(-?\d+(?:\.\d+)?)", file)).group(1)
            vs    = (re.search("vs_(-?\d+(?:\.\d+)?)", file)).group(1)
            vc    = (re.search("vc_(-?\d+(?:\.\d+)?)", file)).group(1)
            vp    = (re.search("vp_(-?\d+(?:\.\d+)?)", file)).group(1)
            params = [float(chips), float(chipc), float(chisc), float(vs), float(vc), float(vp)]
            if params in completed_calculations:
                # print(params)
                same_count += 1
                continue
            else:
                stable_fname = f"chips_{chips}-chipc_{chipc}-chisc_{chisc}-vs_{vs}-vc_{vc}-vp_{vp}.stable_islands.pkl"
                crit_fname   = f"chips_{chips}-chipc_{chipc}-chisc_{chisc}-vs_{vs}-vc_{vc}-vp_{vp}.crits.pkl"

                if stable_fname in l_stable:
                    flist.append(stable_fname)
                    flist.append(crit_fname)
                    files.append(flist) 
                else:
                    print(f"The file {file} has no correpondent file in {args.astable}.")

    print(f"same count = {same_count}")

    vs = []
    vp = []
    vc = []
    chips = []
    chipc = []
    chisc = []
    critpkl     = []
    spkl        = []
    upkl        = []
    binodal_pkl = []
    meshpkl = args.mesh

    for files_list in files:
        chips.append(float((re.search("chips_(-?\d+(?:\.\d+)?)", files_list[0])).group(1)))
        chipc.append(float((re.search("chipc_(-?\d+(?:\.\d+)?)", files_list[0])).group(1)))
        chisc.append(float((re.search("chisc_(-?\d+(?:\.\d+)?)", files_list[0])).group(1)))
        vs.append(float((re.search("vs_(-?\d+(?:\.\d+)?)", files_list[0])).group(1)))
        vc.append(float((re.search("vc_(-?\d+(?:\.\d+)?)", files_list[0])).group(1)))
        vp.append(float((re.search("vp_(-?\d+(?:\.\d+)?)", files_list[0])).group(1)))
        critpkl.append(args.acrits[:-1] +files_list[2])
        spkl.append(args.astable[:-1]+files_list[1])
        upkl.append(args.austble[:-1]+files_list[0]) 
        if args.bpkl is None:
            binodal_pkl.append(f"chips_{chips[-1]}-chipc_{chipc[-1]}-chisc_{chisc[-1]}-vs_{vs[-1]}-vc_{vc[-1]}-vp_{vp[-1]}.binodals.pkl")
        elif args.bpkl[-1] == ".":
                binodal_pkl.append(args.bpkl[:-1]+f"chips_{chips[-1]}-chipc_{chipc[-1]}-chisc_{chisc[-1]}-vs_{vs[-1]}-vc_{vc[-1]}-vp_{vp[-1]}.binodals.pkl")
        else: 
            binodal_pkl.append(args.binodal)

    # print(chisc)

    ncycles = len(vs) // nproc
    vs_list = []
    vp_list = []
    vc_list = []
    chips_list = []
    chipc_list = []
    chisc_list = []
    critpkl_list = []
    spkl_list    = []
    upkl_list    = []
    binodal_list = []
    meshpkl_list = []

    # for i in range(ncycles+1):
    for i in range(nproc):
        vs_list.append(vs[i*ncycles:(i+1)*ncycles])
        vp_list.append(vp[i*ncycles:(i+1)*ncycles])
        vc_list.append(vc[i*ncycles:(i+1)*ncycles])
        chips_list.append(chips[i*ncycles:(i+1)*ncycles])
        chipc_list.append(chipc[i*ncycles:(i+1)*ncycles])
        chisc_list.append(chisc[i*ncycles:(i+1)*ncycles])
        spkl_list.append(spkl[i*ncycles:(i+1)*ncycles])
        upkl_list.append(upkl[i*ncycles:(i+1)*ncycles])
        critpkl_list.append(critpkl[i*ncycles:(i+1)*ncycles])
        binodal_list.append(binodal_pkl[i*ncycles:(i+1)*ncycles])
        meshpkl_list.append([meshpkl]*len(chisc[i*ncycles:(i+1)*ncycles]))

    #print(f"vs_list = {vs_list[0]}")
    # print(f"chisc_list = {chisc_list[0]}")
    print(f"vs_list = {vs_list}", flush=True)
    pool = multiprocessing.Pool(processes=nproc)
    pool.starmap(blibrary.execute_parallel, zip(vs_list, vp_list, vc_list, chisc_list, chips_list, chipc_list, critpkl_list, spkl_list, upkl_list, meshpkl_list, binodal_list))
    print(f"Pool has been closed. This pool had {nproc} threads.", flush=True)

    pool.close()
    pool.join()
    