"""
Validation du théorème d'Ehrenfest.

Règle R4.4 : d⟨R⟩/dt = ⟨P⟩/m et d⟨P⟩/dt = -⟨∇V(R)⟩
Source : [file:1, Chapitre III, § D-1-d]
"""

import numpy as np
from typing import Dict, Callable
from quantum_simulation.core.state import WaveFunctionState
from quantum_simulation.core.operators import PositionOperator, MomentumOperator
from quantum_simulation.utils.numerical import gradient_1d


class EhrenfestValidator:
    """
    Valide théorème d'Ehrenfest : équations classiques pour valeurs moyennes.
    """
    
    def __init__(self, mass: float, hbar: float, tolerance: float = 1e-8):
        """
        Args:
            mass: Masse particule (kg)
            hbar: Constante de Planck réduite (J·s)
            tolerance: Tolérance validation
        """
        self.mass = mass
        self.hbar = hbar
        self.tolerance = tolerance
        
    def compute_mean_values(self, state: WaveFunctionState) -> Dict[str, float]:
        """
        Calcule ⟨X⟩ et ⟨P⟩ pour un état.
        
        Returns:
            Dictionnaire {'mean_x': ..., 'mean_p': ...}
        """
        X_op = PositionOperator(dimension=1)
        P_op = MomentumOperator(hbar=self.hbar, dimension=1)
        
        mean_x = X_op.expectation_value(state)
        mean_p = P_op.expectation_value(state)
        
        return {
            'mean_x': mean_x,
            'mean_p': mean_p
        }
        
    def compute_mean_force(self, state: WaveFunctionState,
                            potential: Callable[[np.ndarray], np.ndarray]) -> float:
        """
        Calcule ⟨-∇V(R)⟩ = ⟨F⟩.
        
        Args:
            state: État quantique
            potential: Fonction potentiel V(x)
            
        Returns:
            Valeur moyenne force (N)
        """
        x_grid = state.spatial_grid
        dx = state.dx
        psi = state.wavefunction
        
        # Valeurs potentiel sur grille
        V_values = potential(x_grid)
        
        # Gradient potentiel (force classique locale)
        grad_V = gradient_1d(V_values, dx, order=2)
        
        # Valeur moyenne : ⟨-∇V⟩ = ∫ ψ*(-∇V)ψ dx
        rho = np.abs(psi)**2
        mean_force = -np.sum(rho * grad_V) * dx
        
        return mean_force
        
    def validate_position_evolution(self, states: list[WaveFunctionState],
                                    times: np.ndarray) -> Dict[str, any]:
        """
        Vérifie d⟨X⟩/dt = ⟨P⟩/m.
        
        Args:
            states: États à différents temps
            times: Temps correspondants
            
        Returns:
            Résultats validation première équation Ehrenfest
        """
        n_times = len(times)
        mean_x = np.zeros(n_times)
        mean_p = np.zeros(n_times)
        
        # Calcul valeurs moyennes
        for i, state in enumerate(states):
            values = self.compute_mean_values(state)
            mean_x[i] = values['mean_x']
            mean_p[i] = values['mean_p']
            
        # Dérivée numérique d⟨X⟩/dt
        dt = times[1] - times[0]  # Supposé uniforme
        dx_dt_numeric = np.gradient(mean_x, dt)
        
        # Côté droit : ⟨P⟩/m
        dx_dt_theory = mean_p / self.mass
        
        # Résidu
        residual = dx_dt_numeric - dx_dt_theory
        max_residual = np.max(np.abs(residual))
        rms_residual = np.sqrt(np.mean(residual**2))
        
        is_satisfied = max_residual < self.tolerance
        
        return {
            'times': times,
            'mean_x': mean_x,
            'mean_p': mean_p,
            'dx_dt_numeric': dx_dt_numeric,
            'dx_dt_theory': dx_dt_theory,
            'residual': residual,
            'max_residual': max_residual,
            'rms_residual': rms_residual,
            'is_satisfied': is_satisfied
        }
        
    def validate_momentum_evolution(self, states: list[WaveFunctionState],
                                    times: np.ndarray,
                                    potential: Callable[[np.ndarray], np.ndarray]) -> Dict[str, any]:
        """
        Vérifie d⟨P⟩/dt = ⟨-∇V(R)⟩.
        
        Returns:
            Résultats validation seconde équation Ehrenfest
        """
        n_times = len(times)
        mean_p = np.zeros(n_times)
        mean_force = np.zeros(n_times)
        
        for i, state in enumerate(states):
            values = self.compute_mean_values(state)
            mean_p[i] = values['mean_p']
            mean_force[i] = self.compute_mean_force(state, potential)
            
        # Dérivée numérique d⟨P⟩/dt
        dt = times[1] - times[0]
        dp_dt_numeric = np.gradient(mean_p, dt)
        
        # Côté droit : ⟨F⟩ = ⟨-∇V⟩
        dp_dt_theory = mean_force
        
        # Résidu
        residual = dp_dt_numeric - dp_dt_theory
        max_residual = np.max(np.abs(residual))
        rms_residual = np.sqrt(np.mean(residual**2))
        
        is_satisfied = max_residual < self.tolerance
        
        return {
            'times': times,
            'mean_p': mean_p,
            'mean_force': mean_force,
            'dp_dt_numeric': dp_dt_numeric,
            'dp_dt_theory': dp_dt_theory,
            'residual': residual,
            'max_residual': max_residual,
            'rms_residual': rms_residual,
            'is_satisfied': is_satisfied
        }
        
    def validate_full_ehrenfest(self, states: list[WaveFunctionState],
                                times: np.ndarray,
                                potential: Callable[[np.ndarray], np.ndarray]) -> Dict[str, any]:
        """
        Validation complète théorème Ehrenfest (position + impulsion).
        
        Returns:
            Rapport complet validation
        """
        position_results = self.validate_position_evolution(states, times)
        momentum_results = self.validate_momentum_evolution(states, times, potential)
        
        overall_valid = (position_results['is_satisfied'] and 
                        momentum_results['is_satisfied'])
        
        return {
            'position_equation': position_results,
            'momentum_equation': momentum_results,
            'overall_valid': overall_valid
        }