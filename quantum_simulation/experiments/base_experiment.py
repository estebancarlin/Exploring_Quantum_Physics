"""
Classe abstraite pour expériences quantiques.

Structure générale : Préparation → Évolution → Mesure → Analyse
"""

from abc import ABC, abstractmethod
from typing import Dict, Any
import time


class Experiment(ABC):
    """
    Classe de base pour toutes expériences de simulation quantique.
    
    Cycle de vie:
        1. prepare_initial_state()
        2. define_hamiltonian()
        3. evolve_state()
        4. perform_measurements()
        5. validate_physics()
        6. analyze_results()
    """
    
    def __init__(self, config: Dict[str, Any]):
        """
        Args:
            config: Dictionnaire configuration (chargé depuis YAML)
        """
        self.config = config
        self.initial_state = None
        self.evolved_states = []
        self.measurement_results = {}
        self.validation_results = {}
        self.execution_time = 0.0
        
    @abstractmethod
    def prepare_initial_state(self):
        """
        Prépare état quantique initial |ψ(t₀)⟩.
        
        Doit vérifier normalisation (Règle R2.1).
        """
        pass
        
    @abstractmethod
    def define_hamiltonian(self):
        """
        Définit hamiltonien H du système.
        
        Doit vérifier hermiticité (Règle R4.5).
        """
        pass
        
    @abstractmethod
    def evolve_state(self):
        """
        Fait évoluer état selon équation Schrödinger (Règle R3.1).
        
        Doit vérifier conservation norme (Règle R5.1).
        """
        pass
        
    @abstractmethod
    def perform_measurements(self):
        """
        Effectue mesures d'observables (Règles R4.1, R4.2).
        """
        pass
        
    @abstractmethod
    def validate_physics(self) -> Dict[str, bool]:
        """
        Valide invariants physiques (Heisenberg, conservation, Ehrenfest).
        
        Returns:
            Dictionnaire {nom_test: réussite}
        """
        pass
        
    def run(self) -> Dict[str, Any]:
        """
        Exécute cycle complet expérience.
        
        Returns:
            Résultats complets (états, mesures, validations, métadonnées)
        """
        start_time = time.time()
        
        print(f"[{self.__class__.__name__}] Démarrage expérience...")
        
        # Étape 1: Préparation
        print("  1/6 Préparation état initial...")
        self.prepare_initial_state()
        
        # Étape 2: Hamiltonien
        print("  2/6 Définition hamiltonien...")
        self.define_hamiltonian()
        
        # Étape 3: Évolution
        print("  3/6 Évolution temporelle...")
        self.evolve_state()
        
        # Étape 4: Mesures
        print("  4/6 Mesures observables...")
        self.perform_measurements()
        
        # Étape 5: Validation
        print("  5/6 Validation physique...")
        self.validation_results = self.validate_physics()
        
        # Étape 6: Analyse (implémentée par sous-classes si nécessaire)
        print("  6/6 Analyse résultats...")
        analysis = self.analyze_results()
        
        self.execution_time = time.time() - start_time
        print(f"  ✓ Expérience terminée en {self.execution_time:.2f}s")
        
        return self._compile_results(analysis)
        
    def analyze_results(self) -> Dict[str, Any]:
        """
        Analyse post-traitement (optionnelle, peut être overridée).
        
        Returns:
            Dictionnaire résultats analyse
        """
        return {}
        
    def _compile_results(self, analysis: Dict[str, Any]) -> Dict[str, Any]:
        """
        Compile tous résultats en structure unique.
        """
        return {
            'experiment_name': self.__class__.__name__,
            'config': self.config,
            'initial_state': self.initial_state,
            'evolved_states': self.evolved_states,
            'measurements': self.measurement_results,
            'validation': self.validation_results,
            'analysis': analysis,
            'execution_time_seconds': self.execution_time,
            'all_validations_passed': all(self.validation_results.values())
        }