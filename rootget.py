import uproot
import hist
import numpy as np
import matplotlib.pyplot as plt


def loadhist(file, treename, filters, bottom, top, bins, label=None, display=None, sample_size = None, loading_step_size = None):
    file = uproot.open(file)
    print(file.classnames())
    tree = file[treename]
    Datapoints = tree.num_entries

    if sample_size == None:
        sample_size = Datapoints
    if loading_step_size == None:
        loading_step_size = int(sample_size/100)

    masshist = hist.Hist(hist.axis.Regular(bins, bottom, top, label=label))

    for baskets in tree.iterate(filter_name = filters, step_size = loading_step_size, entry_stop = sample_size):
        cut = (baskets["nMuon"] == 2)
        phi0 = baskets["Muon_phi", cut, 0]
        phi1 = baskets["Muon_phi", cut, 1]
        eta0 = baskets["Muon_eta", cut, 0]
        eta1 = baskets["Muon_eta", cut, 1]
        pt0 = baskets["Muon_pt", cut, 0]
        pt1 = baskets["Muon_pt", cut, 1]

        mass = np.sqrt(2 * pt0 * pt1 * (np.cosh(eta0 - eta1) - np.cos(phi0 - phi1)))
        masshist.fill(mass)

    if display:
        masshist.plot()
        plt.show()
    return masshist
    
