from quantum_simulation.visualization.quantum_viz_2d import QuantumVisualizer2D
import numpy as np
from Pathlib import Path

def test_density_2d_normalization():
    """Vérifier ∫∫ ρ dxdy = 1."""
    state_2d = create_gaussian_2d(...)
    
    viz = QuantumVisualizer2D()
    fig = viz.plot_density_2d(state_2d)
    
    # Extraction densité depuis figure
    rho = state_2d.probability_density()
    integral = np.trapz(np.trapz(rho, dx), dy)
    
    assert abs(integral - 1.0) < 1e-6
    
def test_animation_2d_creation():
    """Test création animation MP4."""
    states = [create_evolved_state(t) for t in times]
    
    viz = QuantumVisualizer2D()
    output_file = viz.create_animation_2d(states, times, 'test.mp4')
    
    assert Path(output_file).exists()
    assert Path(output_file).stat().st_size > 0