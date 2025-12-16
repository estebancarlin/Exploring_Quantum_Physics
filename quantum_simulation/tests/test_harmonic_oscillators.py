"""
Tests unitaires pour oscillateur harmonique.

Vérifie :
- R6.1 : Spectre Eₙ = ℏω(n + 1/2)
- R6.2 : Algèbre [a, a†] = 1
- R6.3 : Action opérateurs échelle
"""

import sys
from pathlib import Path
# HACK: Ajouter racine projet au path (temporaire)
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import pytest
import numpy as np
from quantum_simulation.systems.harmonic_oscillator import HarmonicOscillator


def test_energy_spectrum():
    """
    Règle R6.1 : Eₙ = ℏω(n + 1/2)
    Source : [file:1, Chapitre V, § B]
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    omega = 1e15
    n_max = 10
    
    ho = HarmonicOscillator(mass, omega, hbar, n_max)
    
    for n in range(n_max + 1):
        E_n = ho.energy_eigenvalue(n)
        E_theory = hbar * omega * (n + 0.5)
        
        assert abs(E_n - E_theory) < 1e-50, \
            f"E_{n} devrait être {E_theory:.2e}, obtenu {E_n:.2e}"


def test_commutator_algebra():
    """
    Règle R6.2 : [a, a†] = 1
    Source : [file:1, Chapitre V, § B]
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    omega = 1e15
    n_max = 50
    
    ho = HarmonicOscillator(mass, omega, hbar, n_max)
    
    # Vérifier algèbre avec tolérance adaptée à troncature
    # Pour n_max=50, erreur bord ~ √(n_max) × eps_machine ≈ 1e-9
    is_valid = ho.validate_algebra(tolerance=1e-8)  # Tolérance relâchée
    assert is_valid, "[a, a†] devrait être identité (modulo troncature)"
    
    # Test alternatif : vérifier commutateur manuellement sur états intérieurs
    a = ho.annihilation_matrix()
    a_dag = ho.creation_matrix()
    commutator = a @ a_dag - a_dag @ a
    identity = np.eye(n_max + 1)
    
    # Vérifier sur états n < n_max (non affectés par troncature)
    for n in range(n_max):
        assert abs(commutator[n, n] - 1.0) < 1e-12, \
            f"[a,a†]_{{n={n}}} devrait être 1"


def test_commutator_algebra_detailed():
    """
    Analyse détaillée erreur commutateur avec troncature.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    omega = 1e15
    
    # Tester différentes troncatures
    for n_max in [10, 20, 50, 100]:
        ho = HarmonicOscillator(mass, omega, hbar, n_max)
        
        a = ho.annihilation_matrix()
        a_dag = ho.creation_matrix()
        commutator = a @ a_dag - a_dag @ a
        identity = np.eye(n_max + 1)
        
        diff = commutator - identity
        
        # Erreur par élément
        max_error = np.max(np.abs(diff))
        diag_error = np.max(np.abs(np.diag(diff)))
        off_diag_error = np.max(np.abs(diff - np.diag(np.diag(diff))))
        
        print(f"\nn_max = {n_max}:")
        print(f"  Erreur max globale  : {max_error:.2e}")
        print(f"  Erreur diagonale    : {diag_error:.2e}")
        print(f"  Erreur hors-diag    : {off_diag_error:.2e}")
        print(f"  Erreur dernier état : {abs(diff[n_max, n_max]):.2e}")
        
        # Commutateur devrait être quasi-parfait sauf dernier état
        assert off_diag_error < 1e-14, "Termes hors-diag devraient être nuls"
        
        # Erreur dernier état attendue : |[a,a†]_{n_max,n_max} - 1| ≈ n_max + 1
        expected_last_error = n_max + 1
        assert abs(abs(diff[n_max, n_max]) - expected_last_error) < 0.1, \
            f"Erreur dernier état devrait être ~{expected_last_error}"
        
        # Erreur intérieure (n < n_max) doit être négligeable
        # diff_interior est DÉJÀ la différence (commutator - identity)
        diff_interior = diff[:n_max, :n_max]
        
        # Vérifier que diff_interior est proche de la matrice identité
        # (c'est-à-dire commutator[:n_max, :n_max] ≈ identity[:n_max, :n_max])
        identity_interior = np.eye(n_max)
        max_error_interior = np.max(np.abs(diff_interior))
        
        assert max_error_interior < 1e-13, \
            f"Erreur intérieure devrait être négligeable pour n_max={n_max}, " \
            f"obtenu {max_error_interior:.2e}"


def test_annihilation_ground_state():
    """
    Règle R6.3 : a|0⟩ = 0 (état fondamental)
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    omega = 1e15
    n_max = 20
    
    ho = HarmonicOscillator(mass, omega, hbar, n_max)
    
    a = ho.annihilation_matrix()
    
    # État fondamental |0⟩ en base {|0⟩, |1⟩, ..., |n_max⟩}
    ground_state = np.zeros(n_max + 1)
    ground_state[0] = 1.0
    
    # a|0⟩
    a_ground = a @ ground_state
    
    # Devrait être vecteur nul
    assert np.linalg.norm(a_ground) < 1e-12, "a|0⟩ devrait être nul"


def test_creation_annihilation_action():
    """
    Règle R6.3 : a|n⟩ = √n|n-1⟩, a†|n⟩ = √(n+1)|n+1⟩
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    omega = 1e15
    n_max = 20
    
    ho = HarmonicOscillator(mass, omega, hbar, n_max)
    
    a = ho.annihilation_matrix()
    a_dag = ho.creation_matrix()
    
    for n in range(1, n_max):
        # État |n⟩
        state_n = np.zeros(n_max + 1)
        state_n[n] = 1.0
        
        # a|n⟩ = √n|n-1⟩
        a_n = a @ state_n
        expected_lower = np.zeros(n_max + 1)
        expected_lower[n-1] = np.sqrt(n)
        
        assert np.allclose(a_n, expected_lower), \
            f"a|{n}⟩ incorrect"
        
        # a†|n⟩ = √(n+1)|n+1⟩
        if n < n_max:
            a_dag_n = a_dag @ state_n
            expected_upper = np.zeros(n_max + 1)
            expected_upper[n+1] = np.sqrt(n+1)
            
            assert np.allclose(a_dag_n, expected_upper), \
                f"a†|{n}⟩ incorrect"


def test_hamiltonian_diagonal():
    """
    Matrice hamiltonien en base propre doit être diagonale.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    omega = 1e15
    n_max = 30
    
    ho = HarmonicOscillator(mass, omega, hbar, n_max)
    
    H = ho.hamiltonian_matrix()
    
    # Vérifier diagonale
    off_diagonal = H - np.diag(np.diag(H))
    assert np.linalg.norm(off_diagonal) < 1e-10, \
        "Hamiltonien devrait être diagonal en base propre"
    
    # Vérifier valeurs diagonales = énergies
    for n in range(n_max + 1):
        E_n = ho.energy_eigenvalue(n)
        assert abs(H[n, n] - E_n) < 1e-48, \
            f"H[{n},{n}] devrait être E_{n}"


def test_number_operator():
    """
    N = a†a doit avoir valeurs propres entières.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    omega = 1e15
    n_max = 15
    
    ho = HarmonicOscillator(mass, omega, hbar, n_max)
    
    a = ho.annihilation_matrix()
    a_dag = ho.creation_matrix()
    N = a_dag @ a
    
    # Vérifier diagonale = n
    for n in range(n_max + 1):
        assert abs(N[n, n] - n) < 1e-10, \
            f"N[{n},{n}] devrait être {n}"
    
    # Vérifier matrice diagonale
    off_diagonal = N - np.diag(np.diag(N))
    assert np.linalg.norm(off_diagonal) < 1e-10, \
        "N devrait être diagonal en base de Fock"


def test_energy_ordering():
    """
    Énergies doivent être croissantes : E₀ < E₁ < E₂ < ...
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    omega = 1e15
    n_max = 100
    
    ho = HarmonicOscillator(mass, omega, hbar, n_max)
    
    energies = [ho.energy_eigenvalue(n) for n in range(n_max + 1)]
    
    for i in range(len(energies) - 1):
        assert energies[i] < energies[i+1], \
            f"E_{i} devrait être < E_{i+1}"


def test_ground_state_energy_positive():
    """
    Énergie fondamentale E₀ = ℏω/2 > 0 (énergie du vide).
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    omega = 1e15
    n_max = 5
    
    ho = HarmonicOscillator(mass, omega, hbar, n_max)
    
    E_0 = ho.energy_eigenvalue(0)
    E_theory = hbar * omega / 2.0
    
    assert E_0 > 0, "Énergie fondamentale doit être positive"
    assert abs(E_0 - E_theory) < 1e-50