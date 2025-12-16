"""
Validation des lois de conservation quantiques.

Règle R5.1 : Conservation de la probabilité totale
Règle R5.2 : Équation de continuité
Source : [file:1, Chapitre III, § D-1-c]
"""

import numpy as np
from typing import Dict, Optional
from quantum_simulation.core.state import WaveFunctionState
from quantum_simulation.utils.numerical import integrate_1d, gradient_1d


class ConservationValidator:
    """
    Valide conservation probabilité et équation de continuité.
    """
    
    def __init__(self, hbar: float, mass: float, tolerance: float = 1e-9):
        """
        Args:
            hbar: Constante de Planck réduite (J·s)
            mass: Masse particule (kg)
            tolerance: Tolérance conservation
        """
        self.hbar = hbar
        self.mass = mass
        self.tolerance = tolerance
        
    def validate_norm_conservation(self, states: list[WaveFunctionState], 
                                    times: np.ndarray) -> Dict[str, any]:
        """
        Vérifie Règle R5.1 : ⟨ψ(t)|ψ(t)⟩ = 1 constant.
        
        Args:
            states: États à différents temps
            times: Temps correspondants
            
        Returns:
            Dictionnaire résultats validation
        """
        norms = np.array([state.norm() for state in states])
        
        # Écarts par rapport à 1
        deviations = np.abs(norms - 1.0)
        max_deviation = np.max(deviations)
        mean_deviation = np.mean(deviations)
        
        is_conserved = max_deviation < self.tolerance
        
        return {
            'times': times,
            'norms': norms,
            'deviations': deviations,
            'max_deviation': max_deviation,
            'mean_deviation': mean_deviation,
            'is_conserved': is_conserved,
            'tolerance': self.tolerance
        }
        
    def compute_probability_current(self, state: WaveFunctionState) -> np.ndarray:
        """
        Calcule courant de probabilité J(x).
        
        Règle R5.2 :  J = (ℏ/2mi)[ψ*∇ψ - ψ∇ψ*]
                        = (1/m)Re(ψ* · (ℏ/i)∇ψ)
        
        Returns:
            Courant J(x) sur grille (m⁻²·s⁻¹ en 3D, m⁻¹·s⁻¹ en 1D)
        """
        psi = state.wavefunction
        dx = state.dx
        
        # Calcul gradients
        grad_psi = gradient_1d(psi, dx, order=2)
        
        # J = (ℏ/2mi)[ψ*∇ψ - ψ∇ψ*]
        # Simplifié : J = (ℏ/m) Im(ψ* ∇ψ)
        psi_conj = np.conj(psi)
        current = (self.hbar / self.mass) * np.imag(psi_conj * grad_psi)
        
        return current
        
    def validate_continuity_equation(self, state_t: WaveFunctionState,
                                    state_t_plus_dt: WaveFunctionState,
                                    dt: float) -> Dict[str, any]:
        """
        Vérifie Règle R5.2 : ∂ρ/∂t + ∇·J = 0
        
        Args:
            state_t: État au temps t
            state_t_plus_dt: État au temps t+dt
            dt: Pas temporel
            
        Returns:
            Résultats validation équation continuité
        """
        # Densités probabilité
        rho_t = state_t.probability_density()
        rho_t_plus_dt = state_t_plus_dt.probability_density()
        
        # Dérivée temporelle (différences finies)
        drho_dt = (rho_t_plus_dt - rho_t) / dt
        
        # Courant moyen entre t et t+dt
        J_t = self.compute_probability_current(state_t)
        J_t_plus_dt = self.compute_probability_current(state_t_plus_dt)
        J_mean = 0.5 * (J_t + J_t_plus_dt)
        
        # Divergence courant
        dx = state_t.dx
        div_J = gradient_1d(J_mean, dx, order=2)
        
        # Équation continuité : résidu = ∂ρ/∂t + ∇·J (devrait être ≈ 0)
        residual = drho_dt + div_J
        
        # Métriques erreur
        max_residual = np.max(np.abs(residual))
        mean_residual = np.mean(np.abs(residual))
        rms_residual = np.sqrt(np.mean(residual**2))
        
        is_satisfied = max_residual < self.tolerance
        
        return {
            'residual': residual,
            'max_residual': max_residual,
            'mean_residual': mean_residual,
            'rms_residual': rms_residual,
            'is_satisfied': is_satisfied,
            'tolerance': self.tolerance,
            'drho_dt': drho_dt,
            'div_J': div_J
        }
        
    def validate_evolution_conservation(self, states: list[WaveFunctionState],
                                        times: np.ndarray) -> Dict[str, any]:
        """
        Validation complète durant évolution : norme + continuité.
        
        Returns:
            Rapport complet conservation
        """
        # Validation norme
        norm_results = self.validate_norm_conservation(states, times)
        
        # Validation continuité (entre temps successifs)
        n_times = len(times)
        continuity_satisfied = np.zeros(n_times - 1, dtype=bool)
        max_residuals = np.zeros(n_times - 1)
        
        for i in range(n_times - 1):
            dt = times[i+1] - times[i]
            cont_result = self.validate_continuity_equation(
                states[i], states[i+1], dt
            )
            continuity_satisfied[i] = cont_result['is_satisfied']
            max_residuals[i] = cont_result['max_residual']
            
        return {
            'norm_conservation': norm_results,
            'continuity_validation': {
                'times_intervals': times[:-1],
                'is_satisfied': continuity_satisfied,
                'max_residuals': max_residuals,
                'all_satisfied': np.all(continuity_satisfied)
            },
            'overall_valid': norm_results['is_conserved'] and np.all(continuity_satisfied)
        }