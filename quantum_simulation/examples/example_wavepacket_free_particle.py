"""
Exemple complet : Évolution paquet d'ondes libre.

Démonstration:
    - Chargement configuration YAML
    - Exécution expérience WavePacketEvolution
    - Validation physique automatique
    - Génération visualisations

Usage:
    python examples/example_wavepacket_free_particle.py
"""

import sys
from pathlib import Path
import yaml

# Ajout chemin projet
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from quantum_simulation.experiments.wavepacket_evolution import WavePacketEvolution
from quantum_simulation.utils.visualization import QuantumVisualizer


def load_configuration(config_path: str = "quantum_simulation/config/parameters.yaml"):
    """
    Charge configuration depuis YAML.
    """
    config_file = project_root / config_path
    
    if not config_file.exists():
        raise FileNotFoundError(f"Fichier config introuvable: {config_file}")
    
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    print(f"Configuration chargée depuis: {config_file}")
    return config


def main():
    """
    Point d'entrée principal.
    """
    print("="*70)
    print(" Simulation Quantique : Évolution Paquet d'Ondes Libre")
    print("="*70)
    print()
    
    # 1. Chargement configuration
    print("[1/4] Chargement configuration...")
    config = load_configuration()
    print()
    
    # 2. Exécution expérience
    print("[2/4] Exécution expérience...")
    experiment = WavePacketEvolution(config)
    results = experiment.run()
    print()
    
    # 3. Affichage résultats validation
    print("[3/4] Résumé validation physique:")
    print("-" * 50)
    for test_name, passed in results['validation'].items():
        if isinstance(passed, bool):
            status = "✓ PASS" if passed else "✗ FAIL"
            print(f"  {test_name:30s} : {status}")
    print("-" * 50)
    
    all_passed = results['all_validations_passed']
    if all_passed:
        print("  ✓ Toutes validations réussies !")
    else:
        print("  ✗ Certaines validations ont échoué")
    print()
    
    # 4. Génération visualisations
    print("[4/4] Génération visualisations...")
    visualizer = QuantumVisualizer(
        output_dir=config.get('visualization', {}).get('output_directory', './results/'),
        dpi=config.get('visualization', {}).get('dpi', 150)
    )
    
    # État initial
    print("  - État initial...")
    visualizer.plot_wavefunction_snapshot(
        state=results['initial_state'],
        time=results['measurements']['times'][0],
        title="État initial - Paquet gaussien",
        save_name="state_initial"
    )
    
    # État final
    print("  - État final...")
    visualizer.plot_wavefunction_snapshot(
        state=results['evolved_states'][-1],
        time=results['measurements']['times'][-1],
        title="État final - Paquet étalé",
        save_name="state_final"
    )
    
    # Évolution observables
    print("  - Évolution observables...")
    visualizer.plot_evolution_observables(
        times=results['measurements']['times'],
        measurement_results=results['measurements'],
        hbar=config['physical_constants']['hbar'],
        save_name="observables_evolution"
    )
    
    # Résumé validation
    print("  - Résumé validation...")
    visualizer.plot_validation_summary(
        validation_results=results['validation'],
        save_name="validation_summary"
    )
    
    # Animation (optionnel)
    if config.get('visualization', {}).get('create_animation', False):
        print("  - Création animation...")
        visualizer.create_evolution_animation(
            states=results['evolved_states'],
            times=results['measurements']['times'],
            interval=200,
            save_name="evolution_animation"
        )
    
    print()
    print("="*70)
    print(f" Simulation terminée en {results['execution_time_seconds']:.2f}s")
    print(f" Résultats sauvegardés dans: {visualizer.output_dir}")
    print("="*70)


if __name__ == "__main__":
    main()