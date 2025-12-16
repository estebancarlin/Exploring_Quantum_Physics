import numpy as np
from quantum_simulation.core.state import EigenStateBasis
from quantum_simulation.core.operators import LadderOperator

class HarmonicOscillator:
    """
    Règles R6.1, R6.2, R6.3
    Source : [file:1, Chapitre V]
    """
    
    def __init__(self, mass: float, omega: float, hbar: float, n_max: int):
        self.mass = mass
        self.omega = omega
        self.hbar = hbar
        self.n_max = n_max
        
    def energy_eigenvalue(self, n: int) -> float:
        """Règle R6.1 : En = ℏω(n + 1/2)"""
        if n < 0 or n > self.n_max:
            raise ValueError(f"n doit être entre 0 et {self.n_max}")
        return self.hbar * self.omega * (n + 0.5)
        
    def hamiltonian_matrix(self) -> np.ndarray:
        """
        Matrice H en base {|0⟩, |1⟩, ..., |n_max⟩}
        Diagonale : En = ℏω(n+1/2)
        """
        H = np.zeros((self.n_max + 1, self.n_max + 1))
        for n in range(self.n_max + 1):
            H[n, n] = self.energy_eigenvalue(n)
        return H
        
    def annihilation_matrix(self) -> np.ndarray:
        """
        Règle R6.3 : a|n⟩ = √n|n-1⟩
        Matrice creuse tri-diagonale inférieure
        """
        a = np.zeros((self.n_max + 1, self.n_max + 1))
        for n in range(1, self.n_max + 1):
            a[n-1, n] = np.sqrt(n)
        return a
        
    def creation_matrix(self) -> np.ndarray:
        """
        Règle R6.3 : a†|n⟩ = √(n+1)|n+1⟩
        Transpose conjuguée de a
        """
        return self.annihilation_matrix().T
        
    def validate_algebra(self, tolerance: float = 1e-10) -> bool:
        """
        Vérifie Règle R6.2 : [a, a†] = I
        """
        a = self.annihilation_matrix()
        a_dag = self.creation_matrix()
        commutator = a @ a_dag - a_dag @ a
        identity = np.eye(self.n_max + 1)
        return np.allclose(commutator, identity, atol=tolerance)
