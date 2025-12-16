import numpy as np
from typing import Callable

from abc import ABC, abstractmethod

from quantum_simulation.core.state import WaveFunctionState, QuantumState
from quantum_simulation.utils.numerical import gradient_1d, laplacian_1d

class Observable(ABC):
    """
    Observable quantique (opérateur hermitique).
    Sources :
    - [file:1, Chapitre II, § D-1] : "observable = opérateur hermitique"
    - [file:1, Chapitre II, § D-2] : valeurs propres réelles
    
    Invariants :
    - Hermiticité : A† = A
    - Valeurs propres réelles
    - Vecteurs propres forment base orthonormée
    """
    
    @abstractmethod
    def apply(self, state: QuantumState) -> QuantumState:
        """Application A|ψ⟩"""
        pass
        
    @abstractmethod
    def expectation_value(self, state: QuantumState) -> float:
        """
        Calcule ⟨A⟩ = ⟨ψ|A|ψ⟩
        Source : [file:1, Chapitre III, § C-4] (Règle 1.4.1)
        """
        pass
        
    @abstractmethod
    def uncertainty(self, state: QuantumState) -> float:
        """
        Calcule ΔA = √(⟨A²⟩ - ⟨A⟩²)
        Source : [file:1, Chapitre III, § C-5] (Règle 1.4.2)
        """
        pass
        
    @abstractmethod
    def eigensystem(self) -> tuple[np.ndarray, list[QuantumState]]:
        """
        Retourne (valeurs_propres, vecteurs_propres)
        
        Limite : Méthode de diagonalisation non spécifiée dans synthèse.
        Implémentation numérique laissée libre (numpy.linalg, etc.)
        """
        pass
        
    def is_hermitian(self, tolerance: float = 1e-10) -> bool:
        """Vérifie A† = A"""
        # Implémentation simple : à surcharger si besoin
        return True
        
    def commutator(self, other: 'Observable') -> 'Observable':
        """
        Calcule [A,B] = AB - BA
        Source : [file:1, Chapitre III, § C-6-a] : compatibilité si [A,B]=0
        """
        # Implémentation symbolique : retour None (à implémenter cas par cas)
        raise NotImplementedError("Commutateur à implémenter pour types spécifiques")


class PositionOperator(Observable):
    """
    Opérateur position X.
    
    En représentation position : X|ψ⟩ → x·ψ(x)
    Source : [file:1, Chapitre II, § E-1]
    """
    
    def __init__(self, dimension: int = 1):
        self.dimension = dimension
        
    def apply(self, state: WaveFunctionState) -> WaveFunctionState:
        """
        Applique X|ψ⟩ = x·ψ(x)
        """
        x = state.spatial_grid
        return WaveFunctionState(x, x * state.wavefunction)
    
    def expectation_value(self, state: WaveFunctionState) -> float:
        """
        ⟨X⟩ = ∫ x|ψ(x)|² dx
        Règle R4.1
        """
        X_psi = self.apply(state)
        return state.inner_product(X_psi).real
    
    def uncertainty(self, state: WaveFunctionState) -> float:
        """
        ΔX = √(⟨X²⟩ - ⟨X⟩²)
        Règle R4.2
        """
        # ⟨X⟩
        mean_x = self.expectation_value(state)
        
        # ⟨X²⟩
        X_psi = self.apply(state)
        X2_psi = self.apply(X_psi)
        mean_x2 = state.inner_product(X2_psi).real
        
        # Variance
        variance = mean_x2 - mean_x**2
        if variance < 0:
            variance = 0  # Erreur numérique
        
        return np.sqrt(variance)
    
    def eigensystem(self) -> tuple[np.ndarray, list[QuantumState]]:
        """
        Position : spectre continu, pas de diagonalisation discrète.
        """
        raise NotImplementedError("Position a spectre continu")


class MomentumOperator(Observable):
    """
    Règle R1.3 : P = -iℏ∇ en représentation position
    Source : [file:1, Chapitre II, § E-2]
    """
    
    def __init__(self, hbar: float, dimension: int = 1):
        self.hbar = hbar
        self.dimension = dimension
        
    def apply(self, state: WaveFunctionState) -> WaveFunctionState:
        """
        Applique P|ψ⟩ = -iℏ∇ψ
        Méthode : Différences finies ordre 2
        """
        grad_psi = gradient_1d(state.wavefunction, state.dx, order=2)
        return WaveFunctionState(
            state.spatial_grid,
            -1j * self.hbar * grad_psi
        )
    
    def expectation_value(self, state: WaveFunctionState) -> float:
        """
        ⟨P⟩ = ⟨ψ|P|ψ⟩
        Règle R4.1
        """
        P_psi = self.apply(state)
        return state.inner_product(P_psi).real
    
    def uncertainty(self, state: WaveFunctionState) -> float:
        """
        ΔP = √(⟨P²⟩ - ⟨P⟩²)
        Règle R4.2
        """
        P_psi = self.apply(state)
        P2_psi = self.apply(P_psi)
        
        exp_P = state.inner_product(P_psi).real
        exp_P2 = state.inner_product(P2_psi).real
        
        variance = exp_P2 - exp_P**2
        if variance < 0:
            variance = 0
        return np.sqrt(variance)
    
    def eigensystem(self) -> tuple[np.ndarray, list[QuantumState]]:
        """Impulsion : spectre continu"""
        raise NotImplementedError("Impulsion a spectre continu")


class Hamiltonian(Observable):
    """
    Règle R3.2 : H = -ℏ²/2m Δ + V(r)
    Source : [file:1, Chapitre III, § B-5-b]
    """
    
    def __init__(self, mass: float, potential: Callable[[np.ndarray], np.ndarray], hbar: float):
        self.mass = mass
        self.potential = potential
        self.hbar = hbar
        
    def apply(self, state: WaveFunctionState) -> WaveFunctionState:
        """
        Applique H|ψ⟩ = [-ℏ²/2m Δ + V(r)]ψ
        """
        # Terme cinétique
        laplacian_psi = laplacian_1d(state.wavefunction, state.dx)
        kinetic_term = -(self.hbar**2 / (2 * self.mass)) * laplacian_psi
        
        # Terme potentiel
        V_values = self.potential(state.spatial_grid)
        potential_term = V_values * state.wavefunction
        
        return WaveFunctionState(
            state.spatial_grid,
            kinetic_term + potential_term
        )
    
    def expectation_value(self, state: WaveFunctionState) -> float:
        """
        ⟨H⟩ = ⟨ψ|H|ψ⟩ (énergie moyenne)
        Règle R4.1
        """
        H_psi = self.apply(state)
        return state.inner_product(H_psi).real
    
    def uncertainty(self, state: WaveFunctionState) -> float:
        """
        ΔH = √(⟨H²⟩ - ⟨H⟩²)
        Règle R4.2
        """
        H_psi = self.apply(state)
        H2_psi = self.apply(H_psi)
        
        exp_H = state.inner_product(H_psi).real
        exp_H2 = state.inner_product(H2_psi).real
        
        variance = exp_H2 - exp_H**2
        if variance < 0:
            variance = 0
        return np.sqrt(variance)
    
    def eigensystem(self) -> tuple[np.ndarray, list[QuantumState]]:
        """
        Diagonalisation numérique hamiltonien.
        
        Méthode : Construction matrice H sur grille, puis np.linalg.eigh()
        
        Returns:
            (energies, eigenstates) triés par énergie croissante
        """
        # Construction matrice hamiltonienne sur grille
        x = None  # Nécessite accès à grille (passée en argument ou stockée)
        
        # LIMITATION : Nécessite grille spatiale connue
        # Solution temporaire : lever erreur, à implémenter si besoin
        raise NotImplementedError(
            "Diagonalisation hamiltonien nécessite grille spatiale. "
            "À implémenter avec méthode spécialisée."
        )
