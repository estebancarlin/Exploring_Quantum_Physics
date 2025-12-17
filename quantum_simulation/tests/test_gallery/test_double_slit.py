from quantum_simulation.experiments.gallery.double_slit_2d import DoubleSlit2D

def test_double_slit_interference():
    """Motif interférence fentes Young."""
    exp = DoubleSlit2D(config)
    results = exp.run()
    
    # Extraction densité écran
    screen_density = results['measurements']['screen_density']
    y_positions = results['measurements']['screen_y']
    
    # Détection pics (scipy.signal.find_peaks)
    peaks, _ = find_peaks(screen_density)
    
    # Interfrange mesuré
    interfrange_measured = np.mean(np.diff(y_positions[peaks]))
    
    # Théorique : Δy = λD/d
    wavelength = hbar / (mass * k0)
    distance_screen = config['screen_distance']
    slit_separation = config['slit_separation']
    interfrange_theory = wavelength * distance_screen / slit_separation
    
    assert abs(interfrange_measured - interfrange_theory) / interfrange_theory < 0.1  # 10%