# Test diagnostic (à exécuter temporairement)

import sys
from pathlib import Path
# HACK: Ajouter racine projet au path (temporaire)
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import numpy as np
from quantum_simulation.systems.potential_systems import InfiniteWell

hbar = 1.054571817e-34
mass = 9.1093837015e-31
width = 1e-9

x = np.linspace(-width, 2*width, 2048)
well = InfiniteWell(width, mass, hbar)

for n in [1, 2, 5, 10]:
    # Version actuelle (sans renormalisation)
    psi = np.zeros_like(x, dtype=complex)
    mask = (x >= 0) & (x <= width)
    norm_factor = np.sqrt(2.0 / width)
    psi[mask] = norm_factor * np.sin(n * np.pi * x[mask] / width)
    
    dx = x[1] - x[0]
    norm = np.sqrt(np.sum(np.abs(psi)**2) * dx)
    print(f"n={n}: norme = {norm:.10f} (écart: {abs(norm-1.0):.2e})")