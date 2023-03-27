import numpy as np
import math
from numba import jit
from scipy.stats import entropy
##########################################################################
#    auxiliary function input qubit spin output it's bloch vector
##########################################################################
# from Cirq
from typing import Iterable, List, Optional, Tuple,Sequence,cast
def validate_indices(num_qubits: int, indices: Sequence[int]) -> None:
    """Validates that the indices have values within range of num_qubits."""
    if any(index < 0 for index in indices):
        raise IndexError(f'Negative index in indices: {indices}')
    if any(index >= num_qubits for index in indices):
        raise IndexError(
            f'Out of range indices, must be less than number of qubits but was {indices}'
        )

def validate_qid_shape(
    state_vector: np.ndarray, qid_shape: Optional[Tuple[int, ...]]
) -> Tuple[int, ...]:
    """Validates the size of the given `state_vector` against the given shape.
    Returns:
        The qid shape.
    Raises:
        ValueError: if the size of `state_vector` does not match that given in
            `qid_shape` or if `qid_state` is not given if `state_vector` does
            not have a dimension that is a power of two.
    """
    size = state_vector.size
    if qid_shape is None:
        qid_shape = (2,) * (size.bit_length() - 1)
    if size != np.prod(qid_shape, dtype=np.int64):
        raise ValueError(
            'state_vector.size ({}) is not a power of two or is not a product '
            'of the qid shape {!r}.'.format(size, qid_shape)
        )
    return qid_shape

def bloch_vector_from_state_vector(
    state_vector: np.ndarray, index: int, qid_shape: Optional[Tuple[int, ...]] = None
) -> np.ndarray:
    """Returns the bloch vector of a qubit.

    Calculates the bloch vector of the qubit at index in the state vector,
    assuming state vector follows the standard Kronecker convention of
    numpy.kron.

    Args:
        state_vector: A sequence representing a state vector in which
            the ordering mapping to qubits follows the standard Kronecker
            convention of numpy.kron (big-endian).
        index: index of qubit who's bloch vector we want to find.
            follows the standard Kronecker convention of numpy.kron.
        qid_shape: specifies the dimensions of the qudits for the input
            `state_vector`.  If not specified, qubits are assumed and the
            `state_vector` must have a dimension a power of two.
            The qudit at `index` must be a qubit.

    Returns:
        A length 3 numpy array representing the qubit's bloch vector.

    Raises:
        ValueError: if the size of `state_vector `is not a power of 2 and the
            shape is not given or if the shape is given and `state_vector` has
            a size that contradicts this shape.
        IndexError: if index is out of range for the number of qubits or qudits
            corresponding to `state_vector`.
    """
    rho = density_matrix_from_state_vector(state_vector, [index], qid_shape=qid_shape)
    v = np.zeros(3, dtype=np.float32)
    v[0] = 2 * np.real(rho[0][1])
    v[1] = 2 * np.imag(rho[1][0])
    v[2] = np.real(rho[0][0] - rho[1][1])


    return v
def von_newman_entropy(rho):
      eigenvalues = np.linalg.eigvalsh(rho)
      return entropy(np.abs(eigenvalues), base=2)
def density_matrix_from_state_vector(
    state_vector: np.ndarray,
    indices: Optional[Iterable[int]] = None,
    qid_shape: Optional[Tuple[int, ...]] = None,
) -> np.ndarray:
    r"""Returns the density matrix of the state vector.

    Calculate the density matrix for the system on the given qubit indices,
    with the qubits not in indices that are present in state vector traced out.
    If indices is None the full density matrix for `state_vector` is returned.
    We assume `state_vector` follows the standard Kronecker convention of
    numpy.kron (big-endian).

    For example:
    state_vector = np.array([1/np.sqrt(2), 1/np.sqrt(2)], dtype=np.complex64)
    indices = None
    gives us

        $$
        \rho = \begin{bmatrix}
                0.5 & 0.5 \\
                0.5 & 0.5
        \end{bmatrix}
        $$

    Args:
        state_vector: A sequence representing a state vector in which
            the ordering mapping to qubits follows the standard Kronecker
            convention of numpy.kron (big-endian).
        indices: list containing indices for qubits that you would like
            to include in the density matrix (i.e.) qubits that WON'T
            be traced out. follows the standard Kronecker convention of
            numpy.kron.
        qid_shape: specifies the dimensions of the qudits for the input
            `state_vector`.  If not specified, qubits are assumed and the
            `state_vector` must have a dimension a power of two.

    Returns:
        A numpy array representing the density matrix.

    Raises:
        ValueError: if the size of `state_vector` is not a power of 2 and the
            shape is not given or if the shape is given and `state_vector`
            has a size that contradicts this shape.
        IndexError: if the indices are out of range for the number of qubits
            corresponding to `state_vector`.
    """
    shape = validate_qid_shape(state_vector, qid_shape)
    n_qubits = len(shape)


    if indices is None:
        return np.outer(state_vector, np.conj(state_vector))


    indices = list(indices)
    validate_indices(n_qubits, indices)


    state_vector = np.asarray(state_vector).reshape(shape)


    sum_inds = np.array(range(n_qubits))
    sum_inds[indices] += n_qubits


    # TODO(#5757): remove type ignore when numpy has proper override signature.
    rho = np.einsum(
        state_vector,
        list(range(n_qubits)),
        np.conj(state_vector),
        cast(List, sum_inds.tolist()),
        indices + cast(List, sum_inds[indices].tolist()),
    )
    new_shape = np.prod([shape[i] for i in indices], dtype=np.int64)


    return rho.reshape((new_shape, new_shape))
def getBlochVector(spin):
    #compute density matrix
    rho = np.outer(spin, np.conj(spin)) #compute density matrix
            
    a = rho[0, 0]
    b = rho[1, 0]
    x = 2.0 * b.real
    y = 2.0 * b.imag
    z = 2.0 * abs(a) - 1.0
    bloch=np.array([x,y,z])        
    return bloch
def pauli4x4Matrixs():    
    p2x = np.array([[0.0, 1.0], [1.0, 0.0]],np.complex128)
    p2y = np.array([[0,-1j],[1j,0]],np.complex128)
    p2z = np.array([[1,0],[0,-1]],np.complex128)

    p4x = np.block([[p2x, np.zeros([2, 2])],                     
                        [np.zeros([2, 2]), p2x]])
    p4y = np.block([[p2y, np.zeros([2, 2])],                     
                        [np.zeros([2, 2]), p2y]])
    p4z = np.block([[p2z, np.zeros([2, 2])],                     
                        [np.zeros([2, 2]), p2z]])
    
    return p4x,p4y,p4z
        
def pauli2x2Matrixs():    
    p2x = np.array([[0.0, 1.0], [1.0, 0.0]],np.complex128)
    p2y = np.array([[0,-1j],[1j,0]],np.complex128)
    p2z = np.array([[1,0],[0,-1]],np.complex128)
    return p2x,p2y,p2z

@jit(nopython=True)
def spinDotBdirac(N: int, DT: float, wf: np.ndarray,B: np.ndarray, spinbloch: np.ndarray):    

    for x in range(N):
        for y in range(N):
            for z in range(N):
                for i in range(4):
                    R = wf[i][x][y][z].real
                    I = wf[i][x][y][z].imag                                            
                    # imag
                    dd = spinbloch.dot(B[x][y][z]) #* prob #(/maxprob)
                    Iinc =  dd * R * DT 
                    #real
                    Rdec = dd * I * DT
                    inc= -Rdec + Iinc *1j
                    #apply B contribution
                    wf[i][x][y][z] += inc    

    
def make2DGaussian(N,L,k,pos,sigma=0.06):
    # p is momentum, k is the wave number that is the# spatial frequency
    # with respect to spatial extent of the simulation
    
    
    p = k * 2.0 * np.pi / L        
    X = L*np.linspace(-0.5, 0.5 - 1.0/N, N)
    Y1, X1 = np.meshgrid(X, X)    
    
    pf = np.exp(1j*p[0]*X1 + 1j*p[1]*Y1)
    gaussian = pf * np.exp(-((X1/L+pos[0]/2.)/sigma)**2/2.0                                    
                  - ((Y1/L-pos[1]/2.)/sigma)**2/2.0 )
    #return gaussian 
    return gaussian/np.sqrt(np.sum(gaussian*np.conj(gaussian)))
    

def getEnergyEigenSpinors(N,L,k,m=1.):
    
    C = 137.036 # Speed of light
    #C=1.0
    mc = m*C
    
    px, py, pz = 2.0*k[0]*np.pi/L, 2.0*k[1]*np.pi/L, 2.0*k[2]*np.pi/L
    
    
    p2 = px**2 + py**2 + pz**2
    p = np.sqrt(p2)
    omega = np.sqrt(mc*mc + p2)

    den1 = p*np.sqrt((mc - omega)**2 + p2)
    den2 = p*np.sqrt((mc + omega)**2 + p2) # from 2d
        
    #zeros = np.zeros(N_DIM*[N], dtype=np.complex128)
    omega = np.sqrt(mc*mc + p2) # Corresponds to E/c

    # Temporary variable used for
    # some denominators in the negative energy solutions
    den1 = p*np.sqrt((mc - omega)**2 + p2)
    # Used for denominators in the positive energy solutions
    den2 = p*np.sqrt((mc + omega)**2 + p2)
    # energy eigenstates in bra form
    # Negative energy solutions
    neg_eig1 = [pz*(mc - omega)/den1, 
                (mc*px - 1.0j*mc*py - (px - 1.0j*py)*omega)/den1,
                p2/den1, 0.]
    neg_eig2 = [(mc*px + 1.0j*mc*py + (-px - 1.0j*py)*omega)/den1,
                -pz*(mc - omega)/den1, 0., p2/den1]

    # Positive energy solutions
    
    pos_eig1 = [pz*(mc + omega)/den2, 
                (mc*px - 1.0j*mc*py + (px - 1.0j*py)*omega)/den2,
                p2/den2, 0.]
    pos_eig2 = [(mc*px + 1.0j*mc*py + (px + 1.0j*py)*omega)/den2,
                -pz*(mc + omega)/den2, 0., p2/den2]
        
    return pos_eig1,pos_eig2, neg_eig1,neg_eig2
def initspinors3D(pos_eig1,pos_eig2,neg_eig1,neg_eig2, spin=None, bPositive=True):    

    if spin is None:
        # if no initial spin is requiered  then spin will point to momentum direction
        spinorsketpos = np.array([pos_eig1[1], 0.0, 0.0, pos_eig1[2]])
        sinfo =" (spin coincides with Kinetic)"
    else:        
        # if initial spin is requiered  make spinors for that specific spin
        # convert spinors to ket form to merge spin to a momentum
        sinfo =" from specific spin" + str(spin) + "\n\t"
        if bPositive:
                            
            sinfo =" positive eigen energies "
            c1 = pos_eig1 @ np.array([spin[0], spin[1], 0.0, 0.0]) # inner product
            c2 = pos_eig2 @ np.array([spin[0], spin[1], 0.0, 0.0])
            spinorsketpos = c1 * np.conj(pos_eig1) + c2 * np.conj(pos_eig2)
        else:
            sinfo =" negative eigen energies "
            c1 = neg_eig1 @ np.array([ 0.0, 0.0,spin[0], spin[1]]) # inner product
            c2 = neg_eig2 @ np.array([ 0.0, 0.0,spin[0], spin[1]])
            spinorsketpos = c1 * np.conj(neg_eig1) + c2 * np.conj(neg_eig2)

    
    spinorsketpos = spinorsketpos / np.linalg.norm(spinorsketpos)        
    print(sinfo)            
    init_spinor = [spinorsketpos[0],spinorsketpos[1],spinorsketpos[2],spinorsketpos[3] ]
    return init_spinor
def makeGaussian2P2D(N,L,pos1,pos2,k1,k2):
    X = L*np.linspace(-0.5, 0.5 - 1.0/N, N)
    Y1, Y2, X1, X2 = np.meshgrid(X, X, X, X)
    k_to_p = (2.0 * np.pi / L)
    px1 = k1[0] * k_to_p
    px2 = k2[0] * k_to_p
    py1 = k1[1] * k_to_p
    py2 = k2[1] * k_to_p
    pf = np.exp(1j*px1*X1 + 1j*py1*Y1 + 1j*px2*X2 + 1j*py2*Y2)
    SIGMA = np.complex128(0.056568)
    xa0=-pos1[0]/2.
    xb0=pos2[0]/2.
    print("xa0",xa0)
    print("xb0",xb0)
    wavefunc = pf * np.exp(-((X1/L+xa0)/SIGMA)**2/2.0
                  - ((X2/L-xb0)/SIGMA)**2/2.0
                  - ((Y1/L)/SIGMA)**2/2.0
                  - ((Y2/L)/SIGMA)**2/2.0
                  )
    wavefunc = wavefunc/np.sqrt(np.sum(wavefunc*np.conj(wavefunc)))
    extent=[Y1[0, 0, 0, 0], Y1[0, -1, 0, 0],
            X1[0, 0, 0, 0], X1[0, 0, -1, 0]]
    print(extent)
    return (wavefunc - np.transpose(wavefunc, (1, 0, 3, 2)))/np.sqrt(2.0), extent
    
def makeGaussian2P2Dold(N,L,pos1,pos2,k1,k2):    

    k_to_p = -(2.0 * np.pi / L)
    px1 = k1[0] * k_to_p
    px2 = k2[0] * k_to_p
    py1 = k1[1] * k_to_p
    py2 = k2[1] * k_to_p
    σ = np.complex128(0.056568)
    
    
    x1, x2,y1,y2 = np.mgrid[ -L/2: L/2:N*1j, -L/2: L/2:N*1j, -L/2: L/2:N*1j, -L/2: L/2:N*1j]

    pf = np.exp(1j*px1*x1 + 1j*py1*y1 + 1j*py2*x2 + 1j*px2*y2)
    
    wf= pf * np.exp( -1/(4* σ**2) * ((x1-pos1[0])**2+(y1-pos1[1])**2+
                                     (x2-pos2[0])**2 +(y2-pos2[1])**2)) / np.sqrt(2*np.pi* σ**2)
    extent=[-1,1,-1,1]

    
    return  wf,extent

    
    
    
