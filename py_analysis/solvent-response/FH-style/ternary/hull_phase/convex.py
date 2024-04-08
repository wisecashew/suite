import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay, delaunay_plot_2d, tsearch, ConvexHull, distance
import scipy.integrate as integrate
import sympy.matrices as smat
import copy
import time

vol_threshold  = 1e-9
conc_threshold = 1e-9
n_mesh         = 251
n_species      = 3

# flory-huggins free energy functional 
def FH_energy(phi, chi, v_molar):
    phi2 = np.maximum(phi, conc_threshold * np.ones_like(phi))
    return np.sum(phi/v_molar * np.log(phi2)) + 1 / 2 * phi.T @ chi @ phi

# map concentration to coordinates
def from_phi_to_xy(phi1_mesh, phi2_mesh):
    return np.array([(2 / np.sqrt(3)) * (1 - phi1_mesh / 2 - phi2_mesh), phi1_mesh]).T

# map coordinates to concentration
def from_xy_to_phi(pt):
    return np.array((pt[:, 1], 1 - (np.sqrt(3) * pt[:, 0] + pt[:, 1]) / 2, (np.sqrt(3) * pt[:, 0] - pt[:, 1]) / 2)).T

# obtain the splits for a particular composition
def calc_phase_split(phi, surface_tot, surf_simplices, cutoff):
    EPS        = 1e-9
    pts_test   = np.append(from_phi_to_xy(phi[0], phi[1]), [0.0])
    simplex_id = -1
    num_simplices = surf_simplices.shape[0]

    for i in range(num_simplices):
        simplex          = surface_tot[surf_simplices[i]] # this is the simplex in question in xy coordinates
        vertices_phi     = from_xy_to_phi(simplex)        # this is the simplex in phi coordinates
        vertices_portion = np.linalg.solve(vertices_phi.T, phi.T)
        if vertices_portion.min() > -EPS:
            simplex_id = i
            break 

    # get the simplex and the A matrix
    simplex      = surface_tot[surf_simplices[simplex_id]]                             # get the simplex again
    A            = 1 * (cutoff > distance.squareform(distance.pdist(simplex[:, :-1]))) # get the A matrix to run the solve

    # first, we distribute the concentration phi to all the vertices (the solution is unique)
    vertices_phi     = from_xy_to_phi(simplex) 
    vertices_portion = np.linalg.solve(vertices_phi, phi.T)

    # then we combine the identical phases
    tot_phase    = 0
    phase_belong = np.ones(n_species) * -1 

    for i in range(n_species):
        if phase_belong[i] < 0:
            phase_belong[i] = tot_phase
            tot_phase += 1
        for j in range(i+1, n_species):
            if A[i, j] > 0:
                phase_belong[j] = phase_belong[i]

    # instantiate containers on how the solution splits 
    phase_composition = np.zeros((n_species, n_species))
    phase_portion     = np.zeros(n_species)
    for i in range(tot_phase):
        for j in range(n_species):
            if phase_belong[j] == i:
                phase_portion[i] += vertices_portion[j]
                phase_composition[i,:] += vertices_phi[j] * vertices_portion[j]
        phase_composition[i, :] /= phase_portion[i]

    ids = np.argsort(phase_composition[:tot_phase,2])
    phase_portion[:tot_phase]        = phase_portion[ids]
    phase_composition[:tot_phase, :] = phase_composition[ids, :]
    phase_portion[tot_phase:]        = np.nan 
    phase_composition[tot_phase:, :] = np.nan 
    return tot_phase, phase_composition, phase_portion, simplex_id


if __name__=="__main__":

    start = time.time()

    # define interaction parameters
    chi = np.array([[0, 2.4, 2.4], [0, 0, 2.4], [0, 0, 0]])
    chi = chi + chi.T

    assert n_species == chi.shape[0]

    # define geometric parameters
    vs = 1
    vp = 1
    vc = 1
    v_molar = np.array([vs, vp, vc])

    # generate the mesh
    phi_list               = np.linspace(0, 1, n_mesh)
    phi_s_mesh, phi_p_mesh = np.meshgrid(phi_list, phi_list, indexing="ij")
    phi_s_mesh = phi_s_mesh.flatten()
    phi_p_mesh = phi_p_mesh.flatten()
    phi_c_mesh = 1.0 - phi_s_mesh - phi_p_mesh

    # curate the compositions
    mask = (phi_c_mesh >= -conc_threshold)
    phi_s_mesh, phi_p_mesh, phi_c_mesh = phi_s_mesh[mask], phi_p_mesh[mask], phi_c_mesh[mask]

    # stack them up
    phi_mesh  = np.vstack((phi_s_mesh, phi_p_mesh, phi_c_mesh)).T 
    mesh_size = phi_mesh.shape[0]

    f_mesh = []
    for i in range(mesh_size):
        free_energy = FH_energy(phi_mesh[i,:], chi, v_molar)
        f_mesh.append(free_energy)
    f_mesh = np.array(f_mesh)
    
    xy_mesh = from_phi_to_xy(phi_mesh[:,0], phi_mesh[:,1])
    surface_raw = np.hstack((xy_mesh, f_mesh.reshape(-1,1)))
    
    # define the caps 
    max_free_energy = 100.0
    cap_free_energy = max_free_energy - 10.0

    surface_lid = np.hstack((xy_mesh, max_free_energy * np.ones((mesh_size,1))))
    
    surface_tot = np.vstack((surface_raw, surface_lid))
    hull        = ConvexHull(surface_tot)

    num_simplices = hull.simplices.shape[0]
    print(f"number of simplices = {num_simplices}")
    simplex_mask  = np.ones(num_simplices, dtype=bool)

    for i in range(num_simplices):
        simplex_pts = surface_tot[hull.simplices[i]]
        if simplex_pts[:, -1].max() >= cap_free_energy:
            simplex_mask[i] = False
        else:
            simplex_pts -= simplex_pts[-1, :]
            simplex_vol  = np.abs(np.linalg.det(simplex_pts[:-1,:-1]))
            if simplex_vol < vol_threshold:
                simplex_mask[i] = False

    surf_simplices = hull.simplices[simplex_mask]
    num_simplices  = surf_simplices.shape[0]

    # get into calculating phases and splits
    nearest_neighbor_distance = (phi_list[1] - phi_list[0]) / np.sqrt(3) * 2
    cutoff                    = nearest_neighbor_distance * 3

    # this is the container for how many phases a certain composition will split into
    num_phases                = np.zeros((num_simplices))

    for i in range(num_simplices):
        simplex_pts   = surface_tot[surf_simplices[i]]
        A             = 1 * (cutoff > distance.squareform(distance.pdist(simplex_pts[:,:-1])))
        D             = np.diag(np.sum(A, axis=1))
        num_phases[i] = smat.Matrix(D-A).eigenvals()[0]

    # plot the phase diagram
    plt.figure(figsize=(10.0, 10.0))
    plt.tripcolor(
        surface_tot[:, 0],
        surface_tot[:, 1],
        surf_simplices,
        facecolors=num_phases,
        edgecolors="k",
    )
    axes = plt.gca()
    axes.set_aspect(1)
    plt.colorbar( label="Number of Phases", orientation="horizontal", fraction=0.046, pad=0.04)
    plt.savefig("FES.png", dpi=1200, bbox_inches="tight")

    print(f"Made FES.png! Moving on to tracing the path...")

    test_size         = 200 
    n_species         = 3
    phi_test          = np.linspace(0, 1, test_size)
    phi_test          = np.vstack((phi_test, (1-phi_test)*2/3, (1-phi_test)/3)).T
    num_phase         = np.zeros(test_size, int)
    phase_composition = np.zeros((test_size, n_species, n_species))
    phase_portion     = np.zeros((test_size, n_species))
    simplex_id        = np.zeros(test_size, dtype=int)

    for i in range(test_size):
        (num_phase[i], phase_composition[i, :, :], phase_portion[i, :], simplex_id[i]) = \
        calc_phase_split(phi_test[i, :], surface_tot, surf_simplices, cutoff)

    # visualize the path in the concentration phase
    pts_test = from_phi_to_xy(phi_test[:,0], phi_test[:,1])

    plt.figure(figsize=(8.0, 8.0))
    plt.tripcolor(
        surface_tot[:, 0],
        surface_tot[:, 1],
        surf_simplices,
        facecolors=num_phases,
        edgecolors="none",
        linewidth=0.01,
    )
    plt.colorbar(label="Number of Phases", orientation="horizontal", fraction=0.046, pad=0.04)
    plt.plot(pts_test[:, 0], pts_test[:, 1], c="r", marker=".", ms=0.5)
    axes = plt.gca()
    axes.set_aspect(1)
    plt.savefig(f"finalFES_.png", dpi=1000, bbox_inches="tight")

    stop = time.time()
    print(f"Time for computation is {stop - start} seconds.")
