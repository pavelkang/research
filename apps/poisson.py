import sys
import cyamites as cy
import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
from scikits.sparse.cholmod import cholesky

def main():
    if(len(sys.argv) > 1):
        filename = sys.argv[1]
    else:
        # TODO Implement a file dialog here
        raise Exception("No file name specified")

    # Load the mesh
    mesh = cy.HalfEdgeMesh(cy.readMesh(filename))
    # 1. Build the matrix

    verts = list(mesh.verts)
    nv = len(verts)
    vd = {verts[i]:i for i in xrange(nv)}

    M_row  = np.array([])
    M_col  = np.array([])
    M_data = np.array([])

    for i in xrange(nv):
        v = verts[i]
        # for M[i,i]
        M_row = np.r_[M_row, [i]]
        M_col = np.r_[M_col, [i]]
        M_data = np.r_[M_data, [1]]
        # for neighbors M[i,j]
        for n in v.adjVerts():
            j = vd[n]
            M_row = np.r_[M_row, [i]]
            M_col = np.r_[M_col, [j]]
            M_data = np.r_[M_data, [1]]
    M = csr_matrix((M_data, (M_row, M_col)), shape=(nv, nv)).todense()

    #print type(mesh.verts)
    # 2. Compute the Cholesky factorization
    # 3. Extract the permutation / reordering of the vertices
    # 4. Visualize the ordering of the vertices

    #print len(mesh.verts)

if __name__ == "__main__":
    main()
