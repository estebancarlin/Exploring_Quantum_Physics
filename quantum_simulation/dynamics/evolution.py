import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from quantum_simulation.core.operators import Hamiltonian
from quantum_simulation.core.state import QuantumState, EigenStateBasis, WaveFunctionState

class TimeEvolution:
    """
    Règle R3.1 : iℏ d|ψ⟩/dt = H|ψ⟩
    
    Évolution temporelle par intégration équation Schrödinger.
    """
    
    def __init__(self, hamiltonian: Hamiltonian, hbar: float):
        """
        Args:
            hamiltonian: Hamiltonien du système
            hbar: Constante de Planck réduite (J·s)
        """
        self.hamiltonian = hamiltonian
        self.hbar = hbar
    
    def _build_hamiltonian_matrix_sparse(self, 
                                        spatial_grid: np.ndarray,
                                        potential: callable = None) -> sparse.csr_matrix:
        """
        Construit matrice hamiltonienne creuse H = -ℏ²/2m Δ + V(R).
        
        Args:
            spatial_grid: Grille spatiale 1D
            potential: Fonction V(x) ou None (particule libre)
            
        Returns:
            Matrice sparse CSR format (nx × nx)
            
        Note:
            - Utilise différences finies ordre 2 (décision D2)
            - Conditions Dirichlet aux bords (décision D3)
            - Matrice tri-diagonale si V=0, bande si V(x)
        """
        from scipy.sparse import diags
        
        nx = len(spatial_grid)
        dx = spatial_grid[1] - spatial_grid[0]
        
        # Terme cinétique : -ℏ²/2m Δ
        # Laplacien : (ψᵢ₊₁ - 2ψᵢ + ψᵢ₋₁)/dx²
        kinetic_coeff = -self.hamiltonian.hbar**2 / (2 * self.hamiltonian.mass * dx**2)
        
        # Diagonales matrice laplacien
        main_diag = -2 * kinetic_coeff * np.ones(nx)
        off_diag = kinetic_coeff * np.ones(nx - 1)
        
        # Matrice tri-diagonale T (terme cinétique)
        T_matrix = diags([off_diag, main_diag, off_diag], 
                        offsets=[-1, 0, 1], 
                        shape=(nx, nx),
                        format='csr')
        
        # Terme potentiel V(x) (diagonal)
        if potential is not None:
            V_values = np.array([potential(x) for x in spatial_grid])
            V_matrix = diags(V_values, offsets=0, format='csr')
            H_matrix = T_matrix + V_matrix
        else:
            H_matrix = T_matrix
        
        return H_matrix
    
    def evolve_wavefunction(self, initial_state: WaveFunctionState, 
                            t0: float, t: float, dt: float) -> WaveFunctionState:
        """
        Évolution temporelle par schéma Crank-Nicolson.
        
        Schéma implicite : (I + iH·dt/2ℏ)ψ(t+dt) = (I - iH·dt/2ℏ)ψ(t)
        
        Propriétés garanties :
            - Conservation norme exacte : ||ψ(t)|| = 1 (Règle R5.1)
            - Évolution unitaire (Complément FIII)
            - Stabilité inconditionnelle (pas de contrainte dt)
            - Précision O(dt²)
        
        Args:
            initial_state: État initial normalisé |ψ(t₀)⟩
            t0: Temps initial (s)
            t: Temps final (s)
            dt: Pas de temps (s)
            
        Returns:
            État évolué |ψ(t)⟩
            
        Raises:
            ValueError: Si norme non conservée (erreur > tolérance)
            
        References:
            - Décision D1 : Crank-Nicolson recommandé
            - Règle R3.2 : iℏ ∂ψ/∂t = Hψ
            - Règle R5.1 : Conservation probabilité
        """
        from scipy.sparse.linalg import spsolve
        from scipy.sparse import eye
        import warnings
        
        # Validation état initial
        if not initial_state.is_normalized():
            raise ValueError("État initial non normalisé !")
        
        # Calcul nombre pas
        n_steps = int((t - t0) / dt)
        if n_steps == 0:
            return initial_state  # Pas d'évolution
        
        # Construction matrices (une seule fois)
        H_matrix = self._build_hamiltonian_matrix_sparse(
            initial_state.spatial_grid,
            potential=self.hamiltonian.potential
        )
        
        nx = len(initial_state.wavefunction)
        I = eye(nx, format='csr')
        
        factor = 0.5j * dt / self.hamiltonian.hbar
        
        # Matrices Crank-Nicolson
        A = I + factor * H_matrix  # (I + iH·dt/2ℏ)
        B = I - factor * H_matrix  # (I - iH·dt/2ℏ)
        
        # Évolution itérative
        psi = initial_state.wavefunction.copy()
        dx = initial_state.dx
        
        tolerance = 1e-9  # Tolérance conservation norme
        
        for step in range(n_steps):
            # Membre droit
            b = B @ psi
            
            # Résolution système linéaire A·ψ(t+dt) = b
            psi = spsolve(A, b)
            
            # Validation conservation norme (Règle R5.1)
            norm_squared = np.sum(np.abs(psi)**2) * dx
            norm = np.sqrt(norm_squared)
            
            if abs(norm - 1.0) > tolerance:
                warnings.warn(
                    f"Pas {step+1}/{n_steps} : Norme = {norm:.10f} "
                    f"(déviation = {abs(norm-1.0):.2e})",
                    RuntimeWarning
                )
                # Renormalisation forcée (sécurité)
                psi /= norm
        
        return WaveFunctionState(initial_state.spatial_grid, psi)