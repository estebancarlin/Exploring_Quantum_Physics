from typing import List, Dict
from quantum_simulation.core.state import WaveFunctionState

class ComparisonEngine:
    """
    Comparaison quantitative entre expériences.
    
    Méthodes:
        - compare_observables(): Écarts relatifs ⟨X⟩, ⟨P⟩, etc.
        - compare_evolutions(): Distances L² entre ψ₁(t) et ψ₂(t)
        - statistical_tests(): Tests χ² pour distributions
        - convergence_analysis(): Études dépendance dt, dx
    """
    
    def compare_observables(self, results_list: List[Dict]) -> ComparisonReport:
        """
        Compare valeurs moyennes observables entre N expériences.
        
        Returns:
            Tableau (N x M) avec M observables + tests significativité
        """
        
    def compare_wavefunctions(self, states1: List[WaveFunctionState],
                            states2: List[WaveFunctionState]) -> Dict:
        """
        Calcule distances entre fonctions d'onde.
        
        Métriques:
            - Distance L² : ∫|ψ₁-ψ₂|² dx
            - Fidélité : |⟨ψ₁|ψ₂⟩|²
            - Entropie relative (si applicable)
        """