"""
Utilitaires numériques pour calculs physiques.

Fonctions :
- Différences finies (gradient, laplacien)
- Intégration numérique
- FFT (pour méthodes alternatives)
"""

import numpy as np
from typing import Literal


def gradient_1d(f: np.ndarray, dx: float, 
                order: Literal[2, 4] = 2,
                boundary: Literal['dirichlet', 'periodic', 'neumann'] = 'dirichlet') -> np.ndarray:
    """
    Calcule gradient df/dx par différences finies.
    
    Args:
        f: Fonction discrète sur grille uniforme
        dx: Pas spatial
        order: Ordre précision (2 ou 4)
        boundary: Conditions limites
            - 'dirichlet': f=0 aux bords (dérivée unilatérale)
            - 'periodic': f périodique
            - 'neumann': df/dx=0 aux bords
    
    Returns:
        Gradient df/dx sur même grille
        
    Notes:
        Ordre 2 : 
            df/dx ≈ (f[i+1] - f[i-1]) / (2dx)    (centré)
            Erreur : O(dx²)
        
        Ordre 4 :
            df/dx ≈ (-f[i+2] + 8f[i+1] - 8f[i-1] + f[i-2]) / (12dx)
            Erreur : O(dx⁴)
    """
    n = len(f)
    grad = np.zeros_like(f)
    
    if order == 2:
        # Différences centrées ordre 2
        if boundary == 'periodic':
            grad[0] = (f[1] - f[-1]) / (2 * dx)
            grad[-1] = (f[0] - f[-2]) / (2 * dx)
            for i in range(1, n-1):
                grad[i] = (f[i+1] - f[i-1]) / (2 * dx)
        else:  # dirichlet ou neumann
            # Bords : différences avant/arrière ordre 1
            grad[0] = (f[1] - f[0]) / dx
            grad[-1] = (f[-1] - f[-2]) / dx
            # Intérieur : centré ordre 2
            for i in range(1, n-1):
                grad[i] = (f[i+1] - f[i-1]) / (2 * dx)
    
    elif order == 4:
        # Différences centrées ordre 4
        if boundary == 'periodic':
            for i in range(n):
                im2 = (i - 2) % n
                im1 = (i - 1) % n
                ip1 = (i + 1) % n
                ip2 = (i + 2) % n
                grad[i] = (-f[ip2] + 8*f[ip1] - 8*f[im1] + f[im2]) / (12 * dx)
        else:
            # Bords : ordre 2 (fallback)
            grad[0] = (f[1] - f[0]) / dx
            grad[1] = (f[2] - f[0]) / (2 * dx)
            grad[-2] = (f[-1] - f[-3]) / (2 * dx)
            grad[-1] = (f[-1] - f[-2]) / dx
            # Intérieur : ordre 4
            for i in range(2, n-2):
                grad[i] = (-f[i+2] + 8*f[i+1] - 8*f[i-1] + f[i-2]) / (12 * dx)
    
    else:
        raise ValueError(f"Ordre {order} non supporté (2 ou 4 uniquement)")
    
    return grad


def laplacian_1d(f: np.ndarray, dx: float,
                 order: Literal[2, 4] = 2,
                 boundary: Literal['dirichlet', 'periodic'] = 'dirichlet') -> np.ndarray:
    """
    Calcule laplacien d²f/dx² par différences finies.
    
    Args:
        f: Fonction discrète
        dx: Pas spatial
        order: Ordre précision
        boundary: Conditions limites
    
    Returns:
        Laplacien d²f/dx²
        
    Notes:
        Ordre 2 :
            d²f/dx² ≈ (f[i+1] - 2f[i] + f[i-1]) / dx²
            Erreur : O(dx²)
        
        Ordre 4 :
            d²f/dx² ≈ (-f[i+2] + 16f[i+1] - 30f[i] + 16f[i-1] - f[i-2]) / (12dx²)
            Erreur : O(dx⁴)
    """
    n = len(f)
    laplacian = np.zeros_like(f)
    
    if order == 2:
        if boundary == 'periodic':
            laplacian[0] = (f[1] - 2*f[0] + f[-1]) / (dx**2)
            laplacian[-1] = (f[0] - 2*f[-1] + f[-2]) / (dx**2)
            for i in range(1, n-1):
                laplacian[i] = (f[i+1] - 2*f[i] + f[i-1]) / (dx**2)
        else:  # dirichlet
            # Bords : f=0 implicite hors grille
            laplacian[0] = (f[1] - 2*f[0]) / (dx**2)
            laplacian[-1] = (-2*f[-1] + f[-2]) / (dx**2)
            # Intérieur
            for i in range(1, n-1):
                laplacian[i] = (f[i+1] - 2*f[i] + f[i-1]) / (dx**2)
    
    elif order == 4:
        if boundary == 'periodic':
            for i in range(n):
                im2 = (i - 2) % n
                im1 = (i - 1) % n
                ip1 = (i + 1) % n
                ip2 = (i + 2) % n
                laplacian[i] = (-f[ip2] + 16*f[ip1] - 30*f[i] + 16*f[im1] - f[im2]) / (12 * dx**2)
        else:
            # Bords : ordre 2 (fallback)
            laplacian[0] = (f[1] - 2*f[0]) / (dx**2)
            laplacian[1] = (f[2] - 2*f[1] + f[0]) / (dx**2)
            laplacian[-2] = (f[-1] - 2*f[-2] + f[-3]) / (dx**2)
            laplacian[-1] = (-2*f[-1] + f[-2]) / (dx**2)
            # Intérieur : ordre 4
            for i in range(2, n-2):
                laplacian[i] = (-f[i+2] + 16*f[i+1] - 30*f[i] + 16*f[i-1] - f[i-2]) / (12 * dx**2)
    
    else:
        raise ValueError(f"Ordre {order} non supporté")
    
    return laplacian


def integrate_1d(f: np.ndarray, dx: float, method: str = 'trapezoid') -> complex:
    """
    Intégration numérique 1D.
    
    Args:
        f: Fonction à intégrer (peut être complexe)
        dx: Pas spatial
        method: 'trapezoid' (trapèzes) ou 'simpson' (Simpson)
    
    Returns:
        Intégrale ∫ f(x) dx
        
    Notes:
        Trapèzes : I ≈ dx * (f[0]/2 + f[1] + ... + f[n-2] + f[n-1]/2)
        Simpson : Nécessite n impair
    """
    if method == 'trapezoid':
        return np.trapz(f, dx=dx)
    
    elif method == 'simpson':
        n = len(f)
        if n % 2 == 0:
            # Nombre pair : utiliser trapèzes sur dernier intervalle
            integral_simpson = integrate_simpson(f[:-1], dx)
            integral_last = (f[-2] + f[-1]) * dx / 2.0
            return integral_simpson + integral_last
        else:
            return integrate_simpson(f, dx)
    
    else:
        raise ValueError(f"Méthode {method} non supportée")


def integrate_simpson(f: np.ndarray, dx: float) -> complex:
    """
    Règle de Simpson 1/3 (nécessite nombre impair de points).
    
    I ≈ dx/3 * (f[0] + 4f[1] + 2f[2] + 4f[3] + ... + 4f[n-2] + f[n-1])
    """
    n = len(f)
    if n < 3:
        return np.trapz(f, dx=dx)
    
    if n % 2 == 0:
        raise ValueError("Simpson nécessite nombre impair de points")
    
    # Coefficients : 1, 4, 2, 4, 2, ..., 4, 1
    coeffs = np.ones(n)
    coeffs[1:-1:2] = 4  # Indices impairs
    coeffs[2:-1:2] = 2  # Indices pairs (sauf extrémités)
    
    return (dx / 3.0) * np.sum(coeffs * f)


def fft_gradient(f: np.ndarray, dx: float) -> np.ndarray:
    """
    Calcule gradient via transformée Fourier (conditions périodiques implicites).
    
    Méthode spectrale : df/dx = F⁻¹[ik F[f]]
    où k sont les nombres d'onde.
    
    Args:
        f: Fonction discrète
        dx: Pas spatial
    
    Returns:
        Gradient df/dx
        
    Note:
        Précision machine (pas d'erreur discrétisation) mais assume périodicité.
    """
    n = len(f)
    
    # Nombres d'onde
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    
    # FFT
    f_hat = np.fft.fft(f)
    
    # Multiplication par ik
    grad_hat = 1j * k * f_hat
    
    # FFT inverse
    grad = np.fft.ifft(grad_hat)
    
    return grad


def fft_laplacian(f: np.ndarray, dx: float) -> np.ndarray:
    """
    Calcule laplacien via FFT.
    
    d²f/dx² = F⁻¹[-k² F[f]]
    
    Args:
        f: Fonction discrète
        dx: Pas spatial
    
    Returns:
        Laplacien d²f/dx²
    """
    n = len(f)
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    
    f_hat = np.fft.fft(f)
    laplacian_hat = -(k**2) * f_hat
    laplacian = np.fft.ifft(laplacian_hat)
    
    return laplacian

def gradient_1d(f: np.ndarray, dx: float, 
                order: Literal[2, 4] = 2,
                boundary: Literal['dirichlet', 'periodic'] = 'dirichlet') -> np.ndarray:
    """
    Calcule gradient ∂f/∂x par différences finies.
    
    Args:
        f: Tableau 1D valeurs fonction
        dx: Pas spatial
        order: Ordre précision (2 ou 4)
        boundary: Conditions limites (décision D3)
        
    Returns:
        Gradient df/dx
        
    Precision:
        - Ordre 2 : O(dx²)
        - Ordre 4 : O(dx⁴)
        
    References:
        - Décision D2 : Ordre 2 par défaut, 4 optionnel
    """
    grad = np.zeros_like(f, dtype=complex)
    
    if order == 2:
        # Schéma centré ordre 2
        if boundary == 'dirichlet':
            # Padding f[0]=f[-1]=0 implicite
            grad[1:-1] = (f[2:] - f[:-2]) / (2*dx)
            grad[0] = (f[1] - 0.0) / (2*dx)  # Bord gauche
            grad[-1] = (0.0 - f[-2]) / (2*dx)  # Bord droit
        elif boundary == 'periodic':
            f_ext = np.concatenate([f[-1:], f, f[:1]])
            grad = (f_ext[2:] - f_ext[:-2]) / (2*dx)
    
    elif order == 4:
        # Schéma centré ordre 4
        if boundary == 'dirichlet':
            # Intérieur
            grad[2:-2] = (-f[4:] + 8*f[3:-1] - 8*f[1:-3] + f[:-4]) / (12*dx)
            
            # Bords : retomber sur ordre 2
            grad[0] = (f[1] - 0.0) / (2*dx)
            grad[1] = (f[2] - f[0]) / (2*dx)
            grad[-2] = (f[-1] - f[-3]) / (2*dx)
            grad[-1] = (0.0 - f[-2]) / (2*dx)
        elif boundary == 'periodic':
            f_ext = np.concatenate([f[-2:], f, f[:2]])
            grad = (-f_ext[4:] + 8*f_ext[3:-1] - 8*f_ext[1:-3] + f_ext[:-4]) / (12*dx)
    
    else:
        raise ValueError(f"Ordre {order} non supporté (2 ou 4)")
    
    return grad