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
    
    def __init__(self, hbar: float, mass: float, tolerance: float = 1e-8):
        """
        Args:
            hbar: Constante de Planck réduite (J·s)
            mass: Masse particule (kg)
            tolerance: Tolérance pour tests conservation
        """
        self.hbar = hbar
        self.mass = mass
        self.tolerance = tolerance
        
    def validate_norm_conservation(self, states, times=None) -> Dict[str, any]:
        """
        Vérifie conservation norme ⟨ψ|ψ⟩ = 1.
        
        Deux modes :
        1. État unique : validate_norm_conservation(state)
        2. Évolution : validate_norm_conservation(states, times)
        
        Args:
            states: État unique OU liste d'états
            times: (optionnel) Temps correspondants si liste états
            
        Returns:
            - Mode simple : {'norm': float, 'is_conserved': bool}
            - Mode évolution : {'norms': array, 'times': array, 'is_conserved': bool, 
                                'max_deviation': float, 'mean_deviation': float}
        """
        # Détection mode : liste ou état unique
        if isinstance(states, list):
            # Mode évolution
            if times is None:
                raise ValueError("Mode évolution nécessite argument 'times'")
            
            return self._validate_norm_evolution(states, times)
        else:
            # Mode état unique
            return self._validate_norm_single(states)
    
    def _validate_norm_single(self, state) -> Dict[str, float]:
        """
        Valide norme d'un seul état.
        
        Returns:
            - 'norm': Norme calculée
            - 'is_conserved': True si |norm - 1| < tolerance
        """
        norm = state.norm()
        is_conserved = abs(norm - 1.0) < self.tolerance
        
        return {
            'norm': norm,
            'is_conserved': is_conserved
        }
    
    def _validate_norm_evolution(self, states: list, times: np.ndarray) -> Dict[str, any]:
        """
        Valide conservation norme durant évolution.
        
        Args:
            states: Liste états quantiques
            times: Temps correspondants
            
        Returns:
            - 'norms': Normes à chaque instant
            - 'times': Temps échantillonnés
            - 'is_conserved': True si toutes normes dans tolérance
            - 'max_deviation': Max|norm - 1|
            - 'mean_deviation': Moyenne|norm - 1|
        """
        if len(states) != len(times):
            raise ValueError(f"Nombre états ({len(states)}) ≠ nombre temps ({len(times)})")
        
        norms = np.array([state.norm() for state in states])
        
        # Déviations par rapport à 1
        deviations = np.abs(norms - 1.0)
        max_deviation = np.max(deviations)
        mean_deviation = np.mean(deviations)
        
        # Conservation : toutes déviations < tolérance
        is_conserved = max_deviation < self.tolerance
        
        return {
            'norms': norms,
            'times': times,
            'is_conserved': is_conserved,
            'max_deviation': max_deviation,
            'mean_deviation': mean_deviation
        }
        
    def compute_probability_current(self, state: WaveFunctionState) -> np.ndarray:
        """
        Calcule courant probabilité J(x) = (ℏ/m) Im(ψ* ∇ψ).
        
        Args:
            state: État quantique
            
        Returns:
            J(x) en kg/s (densité flux probabilité)
            
        Note:
            Formule : J = (ℏ/m) Im(ψ* dψ/dx)
                       = (ℏ/2im) (ψ* dψ/dx - ψ dψ*/dx)
        """
        psi = state.wavefunction
        dx = state.dx
        
        # Gradient de ψ (différences finies centrées)
        dpsi_dx = gradient_1d(psi, dx, order=2)
        
        # Formule : J = (ℏ/m) Im(ψ* ∇ψ)
        #             = (ℏ/m) × (partie imaginaire de ψ* ∇ψ)
        J = (self.hbar / self.mass) * np.imag(np.conj(psi) * dpsi_dx)
        
        return J
        
    def validate_continuity_equation(self, state_t0: WaveFunctionState,
                                    state_t1: WaveFunctionState,
                                    dt: float) -> Dict[str, float]:
        """
        Vérifie équation continuité : ∂ρ/∂t + ∇·J = 0.
        
        Args:
            state_t0: État au temps t
            state_t1: État au temps t + dt
            dt: Pas temporel (s)
            
        Returns:
            - 'max_residual': Max|∂ρ/∂t + ∇·J| (1/s)
            - 'mean_residual': Moyenne|résidu|
            - 'is_satisfied': True si max_residual < tolerance
            
        Note:
            Approximation centrée pour dérivées temporelles :
            ∂ρ/∂t ≈ (ρ(t+dt) - ρ(t)) / dt
            
            Pour divergence :
            ∇·J ≈ dJ/dx (différences finies)
        """
        if not np.allclose(state_t0.spatial_grid, state_t1.spatial_grid):
            raise ValueError("États doivent être sur même grille spatiale")
        
        dx = state_t0.dx
        
        # Densités probabilité
        rho_t0 = state_t0.probability_density()
        rho_t1 = state_t1.probability_density()
        
        # Dérivée temporelle (différence avant)
        drho_dt = (rho_t1 - rho_t0) / dt
        
        # Courant moyen (évaluation à t + dt/2 par moyenne)
        J_t0 = self.compute_probability_current(state_t0)
        J_t1 = self.compute_probability_current(state_t1)
        J_avg = 0.5 * (J_t0 + J_t1)
        
        # Divergence du courant
        div_J = gradient_1d(J_avg, dx, order=2)
        
        # Résidu : ∂ρ/∂t + ∇·J
        residual = drho_dt + div_J
        
        # Métriques
        max_residual = np.max(np.abs(residual))
        mean_residual = np.mean(np.abs(residual))
        
        # Tolérance relative (comparée à échelle typique ∂ρ/∂t)
        drho_dt_scale = np.max(np.abs(drho_dt))
        if drho_dt_scale > 1e-50:
            relative_max_residual = max_residual / drho_dt_scale
            is_satisfied = relative_max_residual < 0.1  # Tolérance 10%
        else:
            # Si ρ stationnaire, résidu doit être petit en absolu
            is_satisfied = max_residual < self.tolerance
        
        return {
            'max_residual': max_residual,
            'mean_residual': mean_residual,
            'relative_max_residual': relative_max_residual if drho_dt_scale > 1e-50 else 0.0,
            'is_satisfied': is_satisfied,
            'drho_dt_scale': drho_dt_scale
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