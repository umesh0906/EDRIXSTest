#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python


# In[115]:


import sys
import numpy as np
import edrixs
from edrixs.fedrixs import ed_fsolver
from mpi4py import MPI


# In[116]:


def get_hopping_coulomb(N_site):
    N = N_site
    norbs = 2*N_site
    U, t = 4.0, -1.0

    hopp = np.zeros((N, N), dtype=np.complex128)
    for i in range(N_site):
        #for j in range(N_site)
        hopp[i, (i-1+N)%N] = hopp[i, (i+1)%N] = t

    hopping = np.zeros((norbs, norbs), dtype=np.complex128)
    hopping[0:norbs:2, 0:norbs:2] = hopp
    hopping[1:norbs:2, 1:norbs:2] = hopp

    umat = np.zeros((norbs, norbs, norbs, norbs), dtype=np.complex128)
    for i in range(N):
        off = i * 2
        umat[off, off + 1, off + 1, off] = U

    edrixs.write_emat(hopping, "hopping_i.in", 1E-10)
    edrixs.write_umat(umat, "coulomb_i.in", 1E-10)


# In[117]:


def get_config():
    config_in = [
        "&control",
        "ed_solver    = 0",
        "num_val_orbs = 2",
        "neval        = 2",
        "nvector      = 2",
        "maxiter      = 500",
        "eigval_tol   = 1E-10",
        "idump        = .true.",
        "&end"
    ]
    f = open('config.in', 'w')
    for line in config_in:
        f.write(line + "\n")
    f.close()


# In[118]:


def get_fock(tot_sz, N_site):
    Sz_list = [1, -1] * N_site
    basis = edrixs.get_fock_basis_by_NSz(2*N_site, N_site, Sz_list)
    print("Total Sz: ", tot_sz)
    for key, val in list(basis.items()):
        if key == tot_sz:
            if len(val) > 0:
                val.sort()
                fname = "fock_i.in"
                f = open(fname, 'w')
                print(len(val), file=f)
                for i in val:
                    print(i, file=f)
                f.close()
            else:
                print("ERROR: Wrong total Sz, check the argument and try again !")
                sys.exit()
            break


# In[125]:


if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    fcomm = comm.py2f()

    if rank == 0:
        N_site =2
        get_hopping_coulomb(N_site)
        get_config()
        tot_sz = 0#int(sys.argv[1])
        get_fock(tot_sz, N_site)
    comm.Barrier()
    print("edrixs >>> Running ED ...")
    ed_fsolver(fcomm, rank, size)
