from quantum_simulation.experiments.base_experiment import Experiment

class TunnelingBarrier(Experiment):
    """
    Effet tunnel à travers barrière potentielle.
    
    Configuration:
        - Barrière rectangulaire : V(x) = V₀ si |x| < L/2, 0 sinon
        - Paquet gaussien incident (E < V₀)
        
    Observables:
        - Coefficient transmission T(E)
        - Coefficient réflexion R(E)
        - Vérification T + R = 1
        
    Scan paramétrique:
        - Variation largeur L (à V₀ fixe)
        - Variation énergie E (à L, V₀ fixes)
        
    Validation:
        - Formule WKB : T ∝ exp[-2∫√(2m(V-E)) dx/ℏ]
    """