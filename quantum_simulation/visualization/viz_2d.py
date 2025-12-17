from typing import List
import numpy as np

class QuantumVisualizer2D:
    """
    Visualisations systèmes 2D (x, y).
    
    Méthodes:
        - plot_density_2d(): Carte densité probabilité ρ(x,y)
        - plot_phase_field(): Carte phase arg(ψ)
        - plot_current_field(): Champ vecteurs J⃗(x,y)
        - plot_potential_overlay(): Potentiel V(x,y) + |ψ|²
        - create_animation_2d(): Animation évolution spatiotemporelle
    """
    
    def plot_density_2d(self, state: WaveFunctionState2D,
                        colormap: str = 'viridis',
                        contour_levels: int = 20,
                       **kwargs) -> Figure:
        """
        Carte 2D densité probabilité.
        
        Args:
            state: État 2D (grille nx × ny)
            colormap: Palette couleurs (perceptuellement uniforme)
            contour_levels: Nombre lignes de niveau
            
        Returns:
            Figure matplotlib avec :
            - Heatmap ρ(x,y) = |ψ(x,y)|²
            - Contours superposés
            - Barre couleur normalisée
            
        Validation:
            - ∫∫ ρ(x,y) dxdy = 1 (affiché légende)
        """
        
    def plot_current_field(self, state: WaveFunctionState2D,
                            hbar: float, mass: float,
                            skip: int = 4) -> Figure:
        """
        Champ vecteurs densité courant J⃗.
        
        Formule (Règle R5.2) :
            J⃗ = (ℏ/2mi)[ψ*∇ψ - ψ∇ψ*]
        
        Args:
            skip: Échantillonnage vecteurs (1 sur skip)
            
        Returns:
            Quiver plot avec :
            - Flèches J⃗ normalisées
            - Couleur = |J⃗| (intensité)
            - Background = ρ(x,y)
        """
        
    def create_animation_2d(self, states: List[WaveFunctionState2D],
                            times: np.ndarray,
                            output_file: str = "evolution_2d.mp4",
                            fps: int = 30) -> str:
        """
        Animation évolution 2D.
        
        Format:
            - MP4 (H.264, compatibilité web)
            - GIF (léger, présentations)
            
        Options:
            - Dual-panel : |ψ|² + arg(ψ)
            - Overlay potentiel V(x,y)
            - Timestamp + observables (⟨X⟩, ⟨Y⟩)
        """