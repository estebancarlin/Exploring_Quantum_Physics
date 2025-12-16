"""
Outils de visualisation pour simulations quantiques.

Génère graphiques pour:
    - Fonctions d'onde ψ(x,t)
    - Densités de probabilité ρ(x,t)
    - Courants de probabilité J(x,t)
    - Évolution observables
    - Résultats validation
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from typing import List, Optional, Dict, Any
from pathlib import Path

from quantum_simulation.core.state import WaveFunctionState


class QuantumVisualizer:
    """
    Générateur de visualisations pour états et évolutions quantiques.
    """
    
    def __init__(self, output_dir: str = "./results/", dpi: int = 150):
        """
        Args:
            output_dir: Dossier sauvegarde figures
            dpi: Résolution figures
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.dpi = dpi
        
        # Style
        plt.style.use('seaborn-v0_8-darkgrid')
        
    def plot_wavefunction_snapshot(self, state: WaveFunctionState, 
                                    time: float,
                                    title: Optional[str] = None,
                                    save_name: Optional[str] = None):
        """
        Affiche |ψ(x)| et Re(ψ), Im(ψ) à un instant donné.
        
        Args:
            state: État quantique
            time: Temps correspondant (pour label)
            title: Titre personnalisé
            save_name: Nom fichier sauvegarde (sans extension)
        """
        x = state.spatial_grid
        psi = state.wavefunction
        
        fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        
        # Subplot 1: Module
        axes[0].plot(x * 1e9, np.abs(psi), 'b-', linewidth=2, label='|ψ(x)|')
        axes[0].fill_between(x * 1e9, 0, np.abs(psi), alpha=0.3)
        axes[0].set_ylabel('|ψ(x)| (m^{-1/2})', fontsize=12)
        axes[0].legend(fontsize=10)
        axes[0].grid(True, alpha=0.3)
        
        # Subplot 2: Partie réelle et imaginaire
        axes[1].plot(x * 1e9, np.real(psi), 'r-', linewidth=2, label='Re(ψ)')
        axes[1].plot(x * 1e9, np.imag(psi), 'g--', linewidth=2, label='Im(ψ)')
        axes[1].axhline(0, color='k', linestyle=':', linewidth=1)
        axes[1].set_xlabel('Position x (nm)', fontsize=12)
        axes[1].set_ylabel('ψ(x) (m^{-1/2})', fontsize=12)
        axes[1].legend(fontsize=10)
        axes[1].grid(True, alpha=0.3)
        
        if title is None:
            title = f'Fonction d\'onde à t = {time:.2e} s'
        fig.suptitle(title, fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        
        if save_name:
            filepath = self.output_dir / f"{save_name}.png"
            plt.savefig(filepath, dpi=self.dpi, bbox_inches='tight')
            print(f"    Figure sauvegardée: {filepath}")
        
        plt.show()
        
    def plot_probability_density(self, state: WaveFunctionState,
                                time: float,
                                save_name: Optional[str] = None):
        """
        Affiche densité de probabilité ρ(x) = |ψ(x)|².
        """
        x = state.spatial_grid
        rho = state.probability_density()
        
        plt.figure(figsize=(10, 6))
        plt.plot(x * 1e9, rho * 1e-9, 'b-', linewidth=2)
        plt.fill_between(x * 1e9, 0, rho * 1e-9, alpha=0.4)
        plt.xlabel('Position x (nm)', fontsize=12)
        plt.ylabel('ρ(x) (nm^{-1})', fontsize=12)
        plt.title(f'Densité de probabilité à t = {time:.2e} s', 
                fontsize=14, fontweight='bold')
        plt.grid(True, alpha=0.3)
        
        # Normalisation check (annotation)
        norm_integral = state.norm()**2
        plt.text(0.98, 0.95, f'∫ρ dx = {norm_integral:.6f}',
                transform=plt.gca().transAxes,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                fontsize=10)
        
        if save_name:
            filepath = self.output_dir / f"{save_name}.png"
            plt.savefig(filepath, dpi=self.dpi, bbox_inches='tight')
        
        plt.show()
        
    def plot_evolution_observables(self, times: np.ndarray,
                                    measurement_results: Dict[str, np.ndarray],
                                    hbar: float,
                                    save_name: Optional[str] = None):
        """
        Évolution temporelle observables: ⟨X⟩, ⟨P⟩, ΔX, ΔP, ⟨E⟩.
        """
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Subplot 1: Position moyenne
        axes[0, 0].plot(times * 1e15, measurement_results['mean_position'] * 1e9, 
                        'b-', linewidth=2, marker='o')
        axes[0, 0].set_xlabel('Temps (fs)', fontsize=11)
        axes[0, 0].set_ylabel('⟨X⟩ (nm)', fontsize=11)
        axes[0, 0].set_title('Position moyenne', fontsize=12, fontweight='bold')
        axes[0, 0].grid(True, alpha=0.3)
        
        # Subplot 2: Impulsion moyenne
        axes[0, 1].plot(times * 1e15, measurement_results['mean_momentum'] * 1e24, 
                        'r-', linewidth=2, marker='s')
        axes[0, 1].set_xlabel('Temps (fs)', fontsize=11)
        axes[0, 1].set_ylabel('⟨P⟩ (×10^{-24} kg·m/s)', fontsize=11)
        axes[0, 1].set_title('Impulsion moyenne', fontsize=12, fontweight='bold')
        axes[0, 1].grid(True, alpha=0.3)
        
        # Subplot 3: Incertitudes
        axes[1, 0].plot(times * 1e15, measurement_results['delta_position'] * 1e9, 
                        'b-', linewidth=2, marker='o', label='ΔX')
        axes[1, 0].set_xlabel('Temps (fs)', fontsize=11)
        axes[1, 0].set_ylabel('ΔX (nm)', fontsize=11)
        axes[1, 0].set_title('Incertitude position', fontsize=12, fontweight='bold')
        axes[1, 0].grid(True, alpha=0.3)
        axes[1, 0].legend()
        
        # Subplot 4: Produit Heisenberg
        heisenberg_bound = hbar / 2.0
        axes[1, 1].plot(times * 1e15, measurement_results['heisenberg_product'] / hbar, 
                        'g-', linewidth=2, marker='^', label='ΔX·ΔP / ℏ')
        axes[1, 1].axhline(0.5, color='k', linestyle='--', linewidth=2, 
                            label='Borne Heisenberg (ℏ/2)')
        axes[1, 1].set_xlabel('Temps (fs)', fontsize=11)
        axes[1, 1].set_ylabel('ΔX·ΔP / ℏ', fontsize=11)
        axes[1, 1].set_title('Produit incertitudes', fontsize=12, fontweight='bold')
        axes[1, 1].grid(True, alpha=0.3)
        axes[1, 1].legend()
        
        plt.tight_layout()
        
        if save_name:
            filepath = self.output_dir / f"{save_name}.png"
            plt.savefig(filepath, dpi=self.dpi, bbox_inches='tight')
        
        plt.show()
        
    def plot_validation_summary(self, validation_results: Dict[str, Any],
                                save_name: Optional[str] = None):
        """
        Résumé visuel validation physique (barres succès/échec).
        """
        tests = []
        results = []
        
        for key, value in validation_results.items():
            if isinstance(value, bool):
                tests.append(key.replace('_', ' ').title())
                results.append(value)
        
        if not tests:
            print("Aucun résultat booléen à afficher")
            return
        
        colors = ['green' if r else 'red' for r in results]
        
        fig, ax = plt.subplots(figsize=(10, 6))
        bars = ax.barh(tests, [1]*len(tests), color=colors, alpha=0.7, edgecolor='black')
        
        # Annotations
        for i, (bar, result) in enumerate(zip(bars, results)):
            label = '✓ PASS' if result else '✗ FAIL'
            ax.text(0.5, i, label, ha='center', va='center', 
                    fontsize=12, fontweight='bold', color='white')
        
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.set_xlabel('Statut validation', fontsize=12)
        ax.set_title('Résultats validation physique', fontsize=14, fontweight='bold')
        ax.invert_yaxis()
        
        plt.tight_layout()
        
        if save_name:
            filepath = self.output_dir / f"{save_name}.png"
            plt.savefig(filepath, dpi=self.dpi, bbox_inches='tight')
        
        plt.show()
        
    def create_evolution_animation(self, states: List[WaveFunctionState],
                                    times: np.ndarray,
                                    interval: int = 100,
                                    save_name: Optional[str] = None):
        """
        Crée animation évolution temporelle ψ(x,t).
        
        Args:
            states: Liste états à différents temps
            times: Temps correspondants
            interval: Délai entre frames (ms)
            save_name: Nom fichier .gif (si sauvegarde)
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        
        x = states[0].spatial_grid * 1e9  # en nm
        line_abs, = ax.plot([], [], 'b-', linewidth=2, label='|ψ(x)|')
        line_re, = ax.plot([], [], 'r--', linewidth=1.5, alpha=0.7, label='Re(ψ)')
        
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(0, 1.2 * np.max([np.max(np.abs(s.wavefunction)) for s in states]))
        ax.set_xlabel('Position x (nm)', fontsize=12)
        ax.set_ylabel('ψ(x)', fontsize=12)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes,
                            verticalalignment='top',
                            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        def init():
            line_abs.set_data([], [])
            line_re.set_data([], [])
            time_text.set_text('')
            return line_abs, line_re, time_text
        
        def animate(i):
            psi = states[i].wavefunction
            line_abs.set_data(x, np.abs(psi))
            line_re.set_data(x, np.real(psi))
            time_text.set_text(f't = {times[i]:.2e} s')
            return line_abs, line_re, time_text
        
        anim = FuncAnimation(fig, animate, init_func=init,
                            frames=len(states), interval=interval, blit=True)
        
        if save_name:
            filepath = self.output_dir / f"{save_name}.gif"
            anim.save(filepath, writer='pillow', fps=1000//interval, dpi=self.dpi//2)
            print(f"    Animation sauvegardée: {filepath}")
        
        plt.show()
        
        return anim