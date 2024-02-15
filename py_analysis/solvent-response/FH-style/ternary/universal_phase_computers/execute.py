import blibrary
import argparse
import warnings
import linecache

parser = argparse.ArgumentParser(description='This will find binodals serially.')
parser.add_argument('--chisc',  metavar='chi_sc',  dest='chi_sc',  type=float,   action='store', help='enter A-C exchange parameter.' )
parser.add_argument('--chips',  metavar='chi_ps',  dest='chi_ps',  type=float,   action='store', help='enter A-B exchange parameter.' )
parser.add_argument('--chipc',  metavar='chi_pc',  dest='chi_pc',  type=float,   action='store', help='enter B-C exchange parameter.' )
parser.add_argument('-vs',      metavar='vs',      dest='vs',      type=float,   action='store', help='specific volume of solvent.')
parser.add_argument('-vc',      metavar='vc',      dest='vc',      type=float,   action='store', help='specific volume of cosolvent.')
parser.add_argument('-vp',      metavar='vp',      dest='vp',      type=float,   action='store', help='specific volume of polymer.')
parser.add_argument('--addr-stable',   metavar='astable',   dest='astable', type=str, action='store', help='enter address of stable islands.')
parser.add_argument('--addr-unstable', metavar='aunstable', dest='austble', type=str, action='store', help='enter address of unstable islands.')
parser.add_argument('--addr-crits',    metavar='acrits',    dest='acrits',  type=str, action='store', help='enter address of critical points.')
parser.add_argument('--mesh',          metavar='mesh',      dest='mesh',    type=str, action='store', help="Name of mesh pkl.")
parser.add_argument('--image-name',    metavar='img',       dest='img',     type=str, action='store', help="Name of image file.", default=None)
parser.add_argument('--no-rtw',        dest='nrtw',      action='store_true',  default=False, help="Don't print out the runtime warning.")
parser.add_argument('--binodal-pkl',   metavar='bpkl',      dest='bpkl',    type=str, action='store', help="Name of binodal pkl.", default=None)
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

    tern_b       = True
    crits_b      = True
    edges_b      = False
    plot_binodal = True
    critpkl      = args.acrits
    spkl         = args.astable
    upkl         = args.austble
    mpkl         = args.mesh
    bpkl         = args.bpkl
    vs = args.vs 
    vc = args.vc 
    vp = args.vp
    chisc = args.chi_sc 
    chips = args.chi_ps
    chipc = args.chi_pc

    blibrary.execute(vs, vp, vc, chisc, chips, chipc, tern_b, crits_b, edges_b, plot_binodal, critpkl, spkl, upkl, mpkl, bpkl, "debug")
