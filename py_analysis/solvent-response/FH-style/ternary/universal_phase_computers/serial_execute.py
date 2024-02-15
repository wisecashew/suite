import os
import re
import blibrary
import argparse
import warnings
import linecache

parser = argparse.ArgumentParser(description='This will find binodals serially.')
parser.add_argument('--addr-stable',   metavar='astable',   dest='astable', type=str, action='store', help='enter address of stable islands.')
parser.add_argument('--addr-unstable', metavar='aunstable', dest='austble', type=str, action='store', help='enter address of unstable islands.')
parser.add_argument('--addr-crits',    metavar='acrits',    dest='acrits',  type=str, action='store', help='enter address of critical points.')
parser.add_argument('--mesh',          metavar='mesh',      dest='mesh',    type=str, action='store', help="Name of mesh pkl.")
parser.add_argument('--image-name',    metavar='img',       dest='img',     type=str, action='store', help="Name of image file.", default=None)
parser.add_argument('--binodal-pkl',   metavar='bpkl',      dest='bpkl',    type=str, action='store', help="Name of binodal pkl.", default=None)
parser.add_argument('--logfile',       metavar='logfile',   dest='logfile', type=str, action='store', help="Name of logfile.", default="binexe.log")
parser.add_argument('--not-ternary',   dest='ntern',   action='store_false', help="The plot to be created is not ternary.", default=True)
parser.add_argument('--no-crits',      dest='ncrit',   action='store_false', help="The plot will not include critical points.", default=True)
parser.add_argument('--no-binodals',   dest='npb',     action='store_false', help="The plot will include binodals.", default=True)
parser.add_argument('--no-rtw',        dest='nrtw',    action='store_true',  help="Dont print out the runtime warning.", default=False)
parser.add_argument('--edges',         dest='edges',   action='store_true' , help="The plot will include edges of the spinodal.", default=False)
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

    files = []

    for file in l_unstable:
        to_del = re.search("TODEL", file)
        if to_del:
            print(f"Ignoring {file}...")
            continue
        else:
            flist = [file]
            chips = (re.search("chips_(-?\d+(?:\.\d+)?)", file)).group(1)
            chipc = (re.search("chipc_(-?\d+(?:\.\d+)?)", file)).group(1)
            chisc = (re.search("chisc_(-?\d+(?:\.\d+)?)", file)).group(1)
            vs    = (re.search("vs_(-?\d+(?:\.\d+)?)", file)).group(1)
            vc    = (re.search("vc_(-?\d+(?:\.\d+)?)", file)).group(1)
            vp    = (re.search("vp_(-?\d+(?:\.\d+)?)", file)).group(1)

            stable_fname = f"chips_{chips}-chipc_{chipc}-chisc_{chisc}-vs_{vs}-vc_{vc}-vp_{vp}.stable_islands.pkl"
            crit_fname   = f"chips_{chips}-chipc_{chipc}-chisc_{chisc}-vs_{vs}-vc_{vc}-vp_{vp}.crits.pkl"

            if stable_fname in l_stable:
                flist.append(stable_fname)
                flist.append(crit_fname)
                files.append(flist) 
            else:
                print(f"The file {file} has no correpondent file in {args.astable}.")

    if args.logfile not in os.listdir("."):

        g = open(args.logfile, 'w')
        for files_list in files:
            chips = float((re.search("chips_(-?\d+(?:\.\d+)?)", files_list[0])).group(1))
            chipc = float((re.search("chipc_(-?\d+(?:\.\d+)?)", files_list[0])).group(1))
            chisc = float((re.search("chisc_(-?\d+(?:\.\d+)?)", files_list[0])).group(1))
            vs    = float((re.search("vs_(-?\d+(?:\.\d+)?)", files_list[0])).group(1))
            vc    = float((re.search("vc_(-?\d+(?:\.\d+)?)", files_list[0])).group(1))
            vp    = float((re.search("vp_(-?\d+(?:\.\d+)?)", files_list[0])).group(1))
            tern_b  = args.ntern
            crit_b  = args.ncrit
            edges_b = args.edges
            plot_b  = args.npb
            critpkl = args.acrits[:-1] +files_list[2]
            spkl    = args.astable[:-1]+files_list[1]
            upkl    = args.austble[:-1]+files_list[0]
            image_name = args.img 
            meshpkl    = args.mesh 
            if args.bpkl is None:
                binodal_pkl = f"chips_{chips}-chipc_{chipc}-chisc_{chisc}-vs_{vs}-vc_{vc}-vp_{vp}.binodals.pkl"
            elif args.bpkl[-1] == ".":
                 binodal_pkl = args.bpkl[:-1]+f"chips_{chips}-chipc_{chipc}-chisc_{chisc}-vs_{vs}-vc_{vc}-vp_{vp}.binodals.pkl"
            else: 
                binodal_pkl = args.binodal
            blibrary.execute(vs, vp, vc, chisc, chips, chipc, \
                            tern_b, crit_b, edges_b, plot_b, \
                                critpkl, spkl, upkl, meshpkl, binodal_pkl, image_name)
            g.write(f"Completed: {files_list[0]}\n")
            g.flush()
        g.close()

    else:
        g = open(args.logfile, 'r')
        pattern = r"Completed: (.+\.pkl)"
        completed = []
        for line in g:
            match = re.search(pattern, line)
            if match:
                completed.append(match.group(1))
        g.close() 

        g = open(args.logfile, 'a')
        for files_list in files:
            if files_list[0] in completed:
                continue 
            else:
                chips = float((re.search("chips_(-?\d+(?:\.\d+)?)", files_list[0])).group(1))
                chipc = float((re.search("chipc_(-?\d+(?:\.\d+)?)", files_list[0])).group(1))
                chisc = float((re.search("chisc_(-?\d+(?:\.\d+)?)", files_list[0])).group(1))
                vs    = float((re.search("vs_(-?\d+(?:\.\d+)?)", files_list[0])).group(1))
                vc    = float((re.search("vc_(-?\d+(?:\.\d+)?)", files_list[0])).group(1))
                vp    = float((re.search("vp_(-?\d+(?:\.\d+)?)", files_list[0])).group(1))
                tern_b  = args.ntern
                crit_b  = args.ncrit
                edges_b = args.edges
                plot_b  = args.npb
                critpkl = args.acrits[:-1] +files_list[2]
                spkl    = args.astable[:-1]+files_list[1]
                upkl    = args.austble[:-1]+files_list[0]
                image_name = args.img 
                meshpkl    = args.mesh 
                if args.bpkl is None:
                    binodal_pkl = f"chips_{chips}-chipc_{chipc}-chisc_{chisc}-vs_{vs}-vc_{vc}-vp_{vp}.binodals.pkl"
                elif args.bpkl[-1] == ".":
                     binodal_pkl = args.bpkl[:-1]+f"chips_{chips}-chipc_{chipc}-chisc_{chisc}-vs_{vs}-vc_{vc}-vp_{vp}.binodals.pkl"
                else: 
                    binodal_pkl = args.binodal
                print(f"stable pkl = {spkl}, unstable pkl = {upkl}, crit pkl = {critpkl}")
                blibrary.execute(vs, vp, vc, chisc, chips, chipc, \
                                tern_b, crit_b, edges_b, plot_b, \
                                    critpkl, spkl, upkl, meshpkl, binodal_pkl, image_name)
                g.write(f"Completed: {files_list[0]}\n")
                g.flush()
        g.close()

