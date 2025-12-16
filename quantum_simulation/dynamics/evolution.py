import numpy as np
from scipy.sparse.linalg import spsolve
from quantum_simulation.core.operators import Hamiltonian
from quantum_simulation.core.state import QuantumState, EigenStateBasis, WaveFunctionState

class TimeEvolution:
    """
    Règle R3.1 : iℏ d|ψ⟩/dt = H|ψ⟩
    """
    
    def evolve_wavefunction(self, initial_state: WaveFunctionState, 
                            t0: float, t: float, dt: float) -> WaveFunctionState:
        """
        Intégration par Crank-Nicolson (décision D1).
        
        Schéma : (I + iH·dt/2ℏ)ψ(t+dt) = (I - iH·dt/2ℏ)ψ(t)
        
        Conserve normalisation (Règle R5.1)
        """
        n_steps = int((t - t0) / dt)
        psi = initial_state.wavefunction.copy()
        
        for step in range(n_steps):
            # Construction matrices (à optimiser : pré-calculer si H constant)
            # Forme : A·ψ_new = B·ψ_old
            # À implémenter : construction matrices creuses via H.apply()
            
            # Résolution système linéaire
            # psi = spsolve(A, B @ psi)
            
            # Vérification conservation norme
            norm_squared = np.sum(np.abs(psi)**2) * initial_state.dx
            if abs(norm_squared - 1.0) > 1e-6:
                raise RuntimeError(f"Norme non conservée : {norm_squared}")
                
        return WaveFunctionState(initial_state.spatial_grid, psi)
