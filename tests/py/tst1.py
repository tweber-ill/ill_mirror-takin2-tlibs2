#
# tlibs2 python interface test
# @author Tobias Weber <tweber@ill.fr>
# @date 9-jun-2020
# @license GPLv3, see 'LICENSE' file
#

import tl2


# get energy transfer from ki and kf
def get_E(ki, kf):
	#E_to_k2 = 2.*co.neutron_mass/hbar_in_meVs**2. / co.elementary_charge*1000. * 1e-20
	E_to_k2 = 0.482596406464	# calculated with scipy, using the formula above

	return (ki**2. - kf**2.) / E_to_k2


datfile = "/home/tw/Projects/repos/skx/exp/data1/elast_1.dat"
dat = tl2.FileInstrBaseD.LoadInstr(datfile)

cnt = dat.GetCountVar()
mon = dat.GetMonVar()
cntcol = dat.GetCol(cnt)
moncol = dat.GetCol(mon)

for point_idx in range(cntcol.size()):
	(h, k, l, ki, kf) = dat.GetScanHKLKiKf(point_idx)
	E = get_E(ki, kf)

	counts = cntcol[point_idx]
	mon_counts = moncol[point_idx]

	print("Q = (%.4f %.4f %.4f), E = %.4f: Monitor: %d, Counts: %d" \
		% (h, k, l, E, mon_counts, counts))
