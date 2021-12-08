# file to calculate the Berry phase and related quantities
# written by Pascal Salzbrenner, pts28@cam.ac.uk
# originally written for CamTB - it is no longer needed there but I put quite a lot of thought into the algorithmic implementation, so on the off-chance that this comes in useful in the future, I'm keeping it here

import numpy as np
from k_points import k_mesh

def berry_phase_calc(k_path, hamiltonian):
    """Function to calculate the Berry phase around the k_path
    :param list k_path: closed path in k-space along which to evaluate the Berry phase
    :param Bloch hamiltonian: an instance of a Hamiltonian of the CamTB Bloch class
    :returns float berry_phase: the Berry phase"""

# currently only works if only a single band is occupied

    vector_products = np.array([], dtype=np.cdouble)

    for i in range(len(k_path)):
        hamiltonian.hamiltonian(k_path[i])
        hamiltonian.solve()

        ground_state_energy = hamiltonian.eigenvalues[0]
        ground_state_vector = hamiltonian.eigenvectors[:,0]

        for j in range(1, len(hamiltonian.eigenvalues)):
            # make sure we take the ground state vector
            if hamiltonian.eigenvalues[j] < ground_state_energy:
                ground_state_energy = hamiltonian.eigenvalues[j]
                ground_state_vector = hamiltonian.eigenvectors[:,j]

        if i != 0:
            vector_products = np.append(vector_products, np.dot(np.conj(past_vector), ground_state_vector))

        past_vector = ground_state_vector

    # discrete Berry phase as defined by Vanderbilt (DOI: 10.1017/9781316662205)
    berry_phase_vector = np.log(vector_products)
    berry_phase = -np.imag(np.sum(berry_phase_vector))

    return berry_phase

def chern_number_calc(hamiltonian, k_points, subspace = None):
    """Function to calculate the Chern number characterising a certain model
    :param Bloch hamiltonian: an instance of a Hamiltonian of the CamTB Bloch class
    :param int k_points: the Brillouin zone is discretised into a k_points * k_points * k_points (depending on the dimension)
    :param list subspace: list of orbitals on which the Chern number is calculated
    mesh"""

    berry_phases = np.array([], dtype=np.double)

    [k_grid, k_disp] = k_mesh(hamiltonian.dimensions, k_points, hamiltonian.k_units,
                    reciprocal_lattice=hamiltonian.reciprocal_lattice)

    # use the k-points in k_grid as the edge of a plate around which to evaluate the berry phase

    for k_point in k_grid:
        # check if we are on the far BZ boundary - in that case, don't create a plate
        if hamiltonian.k_reciprocal:
            if (k_point == 1).any():
                continue
        else:
            continue_flag = False
            for vector in hamiltonian.reciprocal_lattice.values():
                for i in range(len(k_point)):
                    # the lattice vectors will likely have some 0 elements. We are not interested in those
                    if k_point[i] != 0 and np.isclose(k_point[i],vector[i], atol=1e-10):
                        continue_flag = True
                        break
                if continue_flag:
                    break
            if continue_flag:
                continue

        k_path = []
        alter_point = np.copy(k_point)

        # essentially, go along the first vector, then the next, etc; then go in the opposite direction for each vector,
        # to end up in the original place
        for disp in k_disp:
            for i in range(25):
                # 100 seems reasonable for most applications
                k_path.append(np.copy(alter_point))
                alter_point += disp/25

        for disp in k_disp:
            for i in range(25):
                k_path.append(np.copy(alter_point))
                alter_point -= disp/25

        k_path.append(k_point)

        berry_phases = np.append(berry_phases, berry_phase_calc(k_path, hamiltonian))

    chern_number = np.sum(berry_phases)/(2*np.pi)

    if subspace:
        orbitals = ""
        for orbital in subspace:
            orbitals += ("{}_".format(orbital))
        orbitals = orbitals.rstrip(("_"))
    else:
        orbitals = "all"

    with open("chern_number_{}.dat".format(orbitals), "w") as outfile:
        outfile.write("The Chern number in {} orbitals is {}.".format(orbitals, chern_number))

    #TODO not currently numerically stable - only works reliably when k_points is divisible by 10
