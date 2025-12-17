<img src="https://r2cdn.perplexity.ai/pplx-full-logo-primary-dark%402x.png" style="height:64px;margin-right:32px"/>

# Analyse détaillée des décisions techniques D1 à D5

## Document d'aide à la décision technique

**Basé exclusivement sur** :

- Cours de référence : Cohen-Tannoudji, Diu, Laloë - Mécanique Quantique Tome I
- État actuel de l'implémentation : Journal des changements (Janvier 2026)

**Date d'analyse** : 17 décembre 2025

***

## D1 : Schéma d'intégration temporelle pour évolution fonction d'onde

### 1. Contexte physique et numérique

#### Problème physique posé

L'évolution temporelle d'un état quantique |ψ(t)⟩ est régie par l'**équation de Schrödinger dépendante du temps** (Règle R3.1) :

```
iℏ ∂|ψ(t)⟩/∂t = H(t)|ψ(t)⟩
```

En représentation position (Règle R3.2), cette équation devient une équation aux dérivées partielles :

```
iℏ ∂ψ(r,t)/∂t = [-ℏ²/2m Δ + V(r,t)] ψ(r,t)
```

**Source cours** : [file:1, Chapitre I, § B-2 ; Chapitre III, § B-4]

Le cours fournit l'équation formelle mais **aucune méthode de résolution numérique**.

#### Équations concernées du cours

1. **Équation de Schrödinger abstraite** : `iℏ d|ψ⟩/dt = H|ψ⟩`
2. **Forme représentation position** : `iℏ ∂ψ/∂t = Hψ` avec `H = -ℏ²/2m Δ + V(r)`
3. **Solution formelle (cas conservatif)** : `|ψ(t)⟩ = exp(-iH(t-t₀)/ℏ)|ψ(t₀)⟩`

#### Contraintes physiques à respecter

1. **Conservation de la norme** (Règle R5.1) : `⟨ψ(t)|ψ(t)⟩ = constante = 1`
    - **Source** : [file:1, Chapitre III, § D-1-c]
    - Découle de l'hermiticité de H : `d⟨ψ|ψ⟩/dt = (1/iℏ)⟨ψ|[H†-H]|ψ⟩ = 0`
2. **Unitarité de l'évolution** : L'opérateur évolution `U(t,t₀)` doit être unitaire
    - **Source** : [file:1, Complément FIII, § 1-b]
    - Garantit `U†U = UU† = 𝟙`
3. **Équation de continuité** (Règle R5.2) : `∂ρ/∂t + ∇·J = 0`
    - **Source** : [file:1, Chapitre III, § D-1-c]
    - Conséquence locale de la conservation norme

### 2. État actuel de l'implémentation

**Status** : ⚠️ **PARTIEL** (selon journal)

#### Ce qui est déjà implémenté

✅ **Méthodes pour états stationnaires** :

- `TimeEvolution.evolve_eigenstate()` : Décomposition spectrale (Règle R3.3)
    - Formule : `cₙ(t) = cₙ(t₀) exp(-iEₙ(t-t₀)/ℏ)`
- `TimeEvolution.evolve_stationary_state()` : Phase globale (Règle R3.4)
    - Formule : `|ψ(t)⟩ = exp(-iE(t-t₀)/ℏ)|φ⟩`


#### Limites explicites actuelles

❌ **`evolve_wavefunction()` NON fonctionnelle** :

- Structure définie mais **retourne état initial sans évolution**
- **Warning explicite** dans le code (journal, section 1.2)
- Signature définie : `evolve_wavefunction(initial: WaveFunctionState, t0: float, t: float, dt: float) -> WaveFunctionState`


#### Workaround actuel

✅ `WavePacketEvolution` fonctionne via états stationnaires (journal, section 1.2)

- Diagonalisation H préalable
- Évolution par phases `exp(-iEₙt/ℏ)`
- **Limites** : Coûteux si dimension H grande, pas d'évolution continue temps réel


### 3. Options techniques possibles

#### Option 1 : Euler explicite

**Description** : Schéma du premier ordre

```
ψ(t+dt) = ψ(t) - (idt/ℏ)Hψ(t)
```

**Hypothèses physiques implicites** : H constant sur [t, t+dt]

**Avantages** :

- Simplicité extrême d'implémentation
- Coût minimal : 1 application de H par pas

**Inconvénients majeurs** :

- ❌ **Ne conserve PAS la norme** : `||ψ(t+dt)||² ≠ ||ψ(t)||²`
    - Démonstration : `||ψ+δψ||² ≈ 1 + 2Re⟨ψ|δψ⟩ = 1 + 2Re⟨ψ|(-idt/ℏ)Hψ⟩`
    - Comme H hermitique : `⟨ψ|Hψ⟩` réel → `Re(...) = 0` seulement si dt→0
- ❌ **Violation Règle R5.1** (conservation probabilité)
- ❌ **Instabilité numérique** : Croissance exponentielle erreurs
    - Nécessite dt < ℏ/(2Eₘₐₓ) (condition stabilité très restrictive)

**Compatibilité cours** : ❌ Viole contrainte fondamentale conservation norme

**Coût numérique** : O(N) par pas, N = nombre points grille

**Verdict** : **À EXCLURE** - Violation physique majeure

***

#### Option 2 : Crank-Nicolson

**Description** : Schéma implicite ordre 2

```
(𝟙 + iH·dt/2ℏ)ψ(t+dt) = (𝟙 - iH·dt/2ℏ)ψ(t)
```

**Hypothèses physiques implicites** : H évalué au milieu de l'intervalle [t, t+dt]

**Avantages** :

- ✅ **Conserve la norme exactement** : `||ψ(t+dt)|| = ||ψ(t)||`
    - Démonstration : Schéma unitaire, opérateur évolution `U = (1+iHdt/2)⁻¹(1-iHdt/2)`
    - `U†U = [(1-iHdt/2)/(1+iHdt/2)][(1+iHdt/2)/(1-iHdt/2)] = 𝟙`
- ✅ **Respecte Règle R5.1**
- ✅ **Inconditionnellement stable** (pas de contrainte dt)
- ✅ **Précision O(dt²)** (ordre 2)

**Inconvénients** :

- Résolution système linéaire à chaque pas : `Aψ(t+dt) = b`
    - A = matrice `(𝟙 + iH·dt/2ℏ)`
    - b = `(𝟙 - iH·dt/2ℏ)ψ(t)`
- Coût résolution : O(N²) si matrice pleine, O(N) si sparse (H typiquement sparse)

**Compatibilité cours** : ✅ **Parfaitement compatible**

- Conservation norme garantie (Règle R5.1)
- Évolution unitaire (Complément FIII)

**Coût numérique** :

- Construction A : O(N) (matrice tridiagonale si grille 1D)
- Résolution système : O(N) avec solveur sparse (`scipy.sparse.linalg.spsolve`)
- **Total par pas** : O(N)

**Implémentation recommandée** :

```python
from scipy.sparse.linalg import spsolve
# A = I + i*H*dt/(2*hbar) en sparse
# b = (I - i*H*dt/(2*hbar)) @ psi_t
psi_t_plus_dt = spsolve(A, b)
```

**Verdict** : **RECOMMANDÉ PRIORITAIREMENT**

***

#### Option 3 : Split-operator (Trotter-Suzuki)

**Description** : Décomposition exponentielle

```
exp(-iHdt/ℏ) ≈ exp(-iV·dt/2ℏ) · exp(-iT·dt/ℏ) · exp(-iV·dt/2ℏ)
```

où `H = T + V`, `T = P²/2m` (cinétique), `V = V(R)` (potentiel)

**Principe** :

1. Demi-pas potentiel en représentation position : multiplication `ψ → exp(-iV·dt/2ℏ)ψ`
2. FFT → représentation impulsion
3. Pas complet cinétique : multiplication `φ(p) → exp(-ip²dt/2mℏ)φ(p)`
4. FFT⁻¹ → représentation position
5. Demi-pas potentiel : `ψ → exp(-iV·dt/2ℏ)ψ`

**Hypothèses physiques implicites** :

- Approximation Trotter : `exp(A+B) ≈ exp(A)exp(B)` si `[A,B]` petit
- **Conditions périodiques implicites** (FFT)
    - ❌ **Non mentionné dans cours**
    - Impose ψ(x_min) = ψ(x_max) (grille périodique)

**Avantages** :

- ✅ Conservation norme (opérateurs unitaires)
- ✅ Très efficace computationnellement (FFT)
- ✅ Précision O(dt²) (avec décomposition symétrique)
- ✅ Pas de résolution système linéaire

**Inconvénients** :

- ❌ **Conditions périodiques non physiques**
    - Particule "revient" par bord opposé
    - Incompatible avec potentiels confinants (puits)
- ❌ Nécessite implémentation FFT (déjà décision D2 différences finies)
- Erreur Trotter : `||exp(-(iH)dt) - exp(-iTdt)exp(-iVdt)|| = O(dt³·||[T,V]||)`

**Compatibilité cours** : ⚠️ **Partiellement compatible**

- Conservation norme : ✅
- Conditions périodiques : ❌ Hypothèse non justifiée par cours

**Coût numérique** : O(N log N) par pas (FFT domine)

**Verdict** : **Acceptable mais non prioritaire** - Nécessite justification conditions périodiques

***

#### Option 4 : Runge-Kutta 4ème ordre (RK4)

**Description** : Schéma explicite multi-étapes

```
k₁ = -(i/ℏ)H·ψ(t)
k₂ = -(i/ℏ)H·[ψ(t) + dt·k₁/2]
k₃ = -(i/ℏ)H·[ψ(t) + dt·k₂/2]
k₄ = -(i/ℏ)H·[ψ(t) + dt·k₃]
ψ(t+dt) = ψ(t) + (dt/6)(k₁ + 2k₂ + 2k₃ + k₄)
```

**Hypothèses physiques implicites** : H constant sur [t, t+dt]

**Avantages** :

- Précision O(dt⁴) (meilleure que Crank-Nicolson)
- Pas de résolution système linéaire

**Inconvénients majeurs** :

- ❌ **Ne conserve PAS la norme** exactement
    - Erreur accumulation norme : O(dt⁴) par pas, O(T·dt³) sur temps total T
- ❌ **Violation Règle R5.1** (même si petite)
- 4 applications de H par pas (vs 2 pour Crank-Nicolson)

**Compatibilité cours** : ❌ Ne garantit pas conservation norme (contrainte fondamentale)

**Coût numérique** : O(4N) par pas

**Verdict** : **À EXCLURE** - Violation conservation norme, même si faible

***

#### Option 5 : Méthode de Magnus (préservant unitarité)

**Description** : Expansion Magnus pour évolution unitaire exacte

**Problème** : Complexité implémentation très élevée, peu documentée

**Verdict** : **Hors périmètre** - Trop avancé pour contexte actuel

### 4. Évaluation critique

#### Options à exclure

1. **Euler explicite (Option 1)** : ❌ **VIOLATION MAJEURE**
    - Raison : Non-conservation norme → viole Règle R5.1 fondamentale
    - Incompatible avec postulats MQ (interprétation probabiliste)
2. **Runge-Kutta 4 (Option 4)** : ❌ **VIOLATION (même si petite)**
    - Raison : Conservation norme seulement approximative
    - Accumulation erreurs sur temps longs

#### Options acceptables mais non prioritaires

**Split-operator (Option 3)** : ⚠️ **Acceptable avec réserves**

- ✅ Avantages : Performance excellente, conservation exacte
- ❌ Problèmes :
    - Conditions périodiques non justifiées par cours
    - Nécessiterait modification décision D3 (bords grille)
    - Incompatible avec systèmes confinés (puits infini déjà implémenté selon journal)
- **Usage recommandé** : Optimisation future si performance critique, particule libre uniquement


#### Option la plus cohérente

**Crank-Nicolson (Option 2)** : ✅ **OPTION RECOMMANDÉE**

**Justifications** :

1. **Physique** :
    - ✅ Conservation norme **exacte** (Règle R5.1)
    - ✅ Évolution unitaire garantie (Complément FIII)
    - ✅ Équation continuité respectée (Règle R5.2)
2. **Numérique** :
    - ✅ Stabilité inconditionnelle (pas de contrainte dt)
    - ✅ Précision O(dt²) suffisante
    - ✅ Coût O(N) si H sparse (cas général 1D)
3. **Implémentation** :
    - ✅ Compatible architecture actuelle (différences finies D2, bords Dirichlet D3)
    - ✅ Bibliothèques Python robustes (`scipy.sparse`)
    - ✅ Validation simple : vérifier `||ψ(t)|| = 1` à chaque pas
4. **Cohérence cours** :
    - ✅ Respecte tous invariants physiques identifiés
    - ✅ Pas d'hypothèse externe au cours

### 5. Recommandation

#### Décision recommandée

**Implémenter Crank-Nicolson** dans `dynamics/evolution.py`, méthode `evolve_wavefunction()`

#### Justification synthétique

**Physique** : Seule méthode garantissant conservation norme exacte (contrainte fondamentale Règle R5.1) sans hypothèses non justifiées par cours.

**Numérique** : Stabilité inconditionnelle + coût O(N) raisonnable avec solveurs sparse.

**Pratique** : Compatible avec toutes décisions existantes (D2, D3), validation triviale.

#### Pseudo-algorithme

```python
def evolve_wavefunction(self, initial_state: WaveFunctionState, 
                       t0: float, t: float, dt: float) -> WaveFunctionState:
    """
    Évolution Crank-Nicolson : (I + iH·dt/2ℏ)ψ(t+dt) = (I - iH·dt/2ℏ)ψ(t)
    Respecte Règle R5.1 : conservation norme exacte.
    """
    n_steps = int((t - t0) / dt)
    psi = initial_state.wavefunction.copy()
    
    # Construire matrices sparse (une seule fois)
    H_matrix = self._build_hamiltonian_matrix(initial_state.spatial_grid)
    factor = 1j * dt / (2 * self.hbar)
    
    A = sparse.eye(len(psi)) + factor * H_matrix  # (I + iH·dt/2ℏ)
    B_op = sparse.eye(len(psi)) - factor * H_matrix  # (I - iH·dt/2ℏ)
    
    for step in range(n_steps):
        b = B_op @ psi
        psi = spsolve(A, b)
        
        # Validation conservation norme (Règle R5.1)
        norm = np.sqrt(np.sum(np.abs(psi)**2) * dx)
        assert abs(norm - 1.0) < tolerance, "Norme non conservée!"
    
    return WaveFunctionState(initial_state.spatial_grid, psi)
```


#### Impact sur autres décisions

**Couplages directs** :

1. **D2 (gradient/laplacien)** : ✅ Compatible
    - H matrice construite avec différences finies ordre 2 (déjà implémenté)
    - Pas de modification nécessaire
2. **D3 (bords grille)** : ✅ Compatible
    - Conditions Dirichlet (ψ=0 aux bords) naturellement gérées dans matrice H
    - Padding déjà implémenté dans `utils/numerical.py` (journal)
3. **D4 (oscillateur)** : ⚠️ Attention
    - Si oscillateur en représentation position (actuellement base abstraite)
    - Nécessiterait polynômes Hermite (Extension future journal, section 5.1)
4. **D5 (mesures)** : ✅ Indépendant
    - Évolution et mesure découplées (architecture)

**Dépendances indirectes** :

- Tests validation (section 6.2 journal) : ⚠️ Adapter tests conservation norme
    - Actuellement ~95% précision équation continuité
    - Crank-Nicolson devrait améliorer à ~100%


#### Priorité et effort

**Priorité** : 🔴 **HAUTE** (bloqueur critique selon journal, section 1.2)

**Effort estimé** :

- Construction matrices H sparse : 2-3h
- Implémentation schéma Crank-Nicolson : 1-2h
- Tests validation (conservation norme, Ehrenfest) : 2-3h
- **Total** : ~6-8h développement

**Tests requis** :

1. Conservation norme : `||ψ(t)|| = 1` pour tout t
2. Équation continuité : `∂ρ/∂t + ∇·J = 0` (améliorer 95%→100%)
3. Ehrenfest : `d⟨X⟩/dt = ⟨P⟩/m` sur paquet gaussien libre
4. Convergence : vérifier erreur O(dt²) en comparant dt, dt/2, dt/4

***

## D2 : Calcul gradient/laplacien (opérateur impulsion)

### 1. Contexte physique et numérique

#### Problème physique posé

En représentation position, l'**opérateur impulsion** agit comme dérivée (Règle 1.7.1) :

```
P = -iℏ∇
```

**Source cours** : [file:1, Chapitre II, § E-2]

Application sur fonction d'onde : `Pψ(r) = -iℏ∇ψ(r)`

Le **hamiltonien** contient terme cinétique (Règle 1.7.2) :

```
H = P²/2m + V(R) = -ℏ²/2m Δ + V(r)
```

**Source cours** : [file:1, Chapitre III, § B-5-b]

où `Δ = ∇² = ∂²/∂x² + ∂²/∂y² + ∂²/∂z²` (laplacien)

Le cours donne **formules continues uniquement**, aucune discrétisation.

#### Contraintes physiques à respecter

1. **Relations de commutation canoniques** (Règle R1.3) :

```
[Rᵢ, Pⱼ] = iℏδᵢⱼ
```

**Source** : [file:1, Chapitre III, § B-5-a]
2. **Hermiticité de P** (Règle R4.5) : `P† = P`
    - **Source** : [file:1, Chapitre II, § D-1]
    - Garantit valeurs propres réelles
3. **Conservation probabilité** : Dérivées doivent respecter équation continuité
    - Courant probabilité : `J = (ℏ/2mi)[ψ*∇ψ - ψ∇ψ*]` (Règle R5.2)

### 2. État actuel de l'implémentation

**Status** : ✅ **RÉSOLU** (selon journal, section 1.1)

#### Ce qui est implémenté

✅ **Différences finies ordre 2 adoptées** :

**Gradient** (journal : `utils/numerical.py` complet) :

```
∂ψ/∂x ≈ (ψᵢ₊₁ - ψᵢ₋₁)/(2dx)  # Différences finies centrées
```

**Laplacien** :

```
∂²ψ/∂x² ≈ (ψᵢ₊₁ - 2ψᵢ + ψᵢ₋₁)/dx²
```

✅ **Validation** (journal, section 1.1) :

- Hermiticité vérifiée et testée
- Commutateurs `[X,P] = iℏ` fonctionnels
- Tests sur cas connus (onde plane, gaussienne)


#### Limites identifiées

⚠️ **Limite N7** (journal, section 4.2) :

- Erreur O(dx²) visible pour petits σₓ < dx
- États localisés nécessitent grilles fines
- **Solution future** : Ordre 4 optionnel


### 3. Options techniques possibles

#### Option 1 : Différences finies ordre 2 (IMPLÉMENTÉ)

**Description** :

**Gradient centré** :

```
(∇ψ)ᵢ = (ψᵢ₊₁ - ψᵢ₋₁)/(2dx)
```

**Laplacien** :

```
(Δψ)ᵢ = (ψᵢ₊₁ - 2ψᵢ + ψᵢ₋₁)/dx²
```

**Hypothèses physiques** :

- ψ suffisamment lisse (C² au moins)
- dx << longueur d'onde caractéristique

**Avantages** :

- ✅ Simplicité implémentation
- ✅ Hermiticité préservée (schéma centré)
- ✅ Coût O(N) faible
- ✅ Précision O(dx²) suffisante pour la plupart des cas

**Inconvénients** :

- ⚠️ Erreur O(dx²) peut être insuffisante pour états très localisés
- Nécessite grilles fines si σₓ petit

**Compatibilité cours** : ✅ **Parfaite** - Approximation numérique standard de `∇` continu

**Coût numérique** : O(N) par application

**Tests effectués** (selon journal) :

- ✅ Commutateurs [X,P] = iℏ vérifiés
- ✅ Hermiticité P† = P vérifiée
- ✅ Conservation probabilité ~95% (équation continuité)

***

#### Option 2 : Différences finies ordre 4

**Description** :

**Gradient** :

```
(∇ψ)ᵢ = [-ψᵢ₊₂ + 8ψᵢ₊₁ - 8ψᵢ₋₁ + ψᵢ₋₂]/(12dx)
```

**Laplacien** :

```
(Δψ)ᵢ = [-ψᵢ₊₂ + 16ψᵢ₊₁ - 30ψᵢ + 16ψᵢ₋₁ - ψᵢ₋₂]/(12dx²)
```

**Avantages** :

- ✅ Précision O(dx⁴) (amélioration significative)
- ✅ Réduit erreurs pour états localisés (Limite N7)

**Inconvénients** :

- Plus complexe aux bords (nécessite ψᵢ₊₂, ψᵢ₋₂)
- Coût légèrement supérieur
- Hermiticité plus délicate à garantir

**Compatibilité cours** : ✅ Toujours approximation de `∇` continu

**Coût numérique** : O(N) par application (constant légèrement plus élevé)

**Verdict** : **Extension future** - Configurable si besoin précision accrue

***

#### Option 3 : FFT (représentation impulsion)

**Description** : Utiliser transformée Fourier

**Principe** :

1. FFT : ψ(x) → φ(k)
2. Multiplication : `(∇ψ)(k) = ik·φ(k)`, `(Δψ)(k) = -k²·φ(k)`
3. FFT⁻¹ : → représentation position

**Avantages** :

- ✅ Précision spectrale (erreur machine si ψ lisse)
- ✅ Très efficace pour grandes grilles

**Inconvénients** :

- ❌ **Conditions périodiques implicites** ψ(x_min) = ψ(x_max)
- ❌ Incompatible avec décision D3 (bords Dirichlet)
- ❌ Surcoût FFT : O(N log N) vs O(N)

**Compatibilité cours** : ⚠️ Conditions périodiques non justifiées

**Coût numérique** : O(N log N) par application

**Verdict** : **Réservé split-operator** (si implémenté Option 3 de D1)

### 4. Évaluation critique

#### Options à exclure

**Aucune** - Toutes options valides physiquement

#### Options acceptables

1. **Ordre 2 (Option 1)** : ✅ **IMPLÉMENTÉ** - Équilibre optimal simplicité/précision
2. **Ordre 4 (Option 2)** : ✅ **Extension future** - Si précision critique
3. **FFT (Option 3)** : ⚠️ **Réservé cas spéciaux** - Nécessite conditions périodiques

#### Option actuelle cohérente

**Différences finies ordre 2** : ✅ **CONFIRMÉ COMME OPTIMAL**

**Raisons** :

- Implémentation robuste déjà validée (journal, tests complets)
- Précision O(dx²) suffisante pour cas d'usage actuels
- Hermiticité garantie (tests passent)
- Coût O(N) minimal


### 5. Recommandation

#### Décision recommandée

**Conserver différences finies ordre 2** avec **option ordre 4 configurable** future

#### Justification synthétique

**Physique** : Hermiticité vérifiée, commutateurs respectés (Règles R1.3, R4.5)

**Numérique** : Équilibre optimal précision O(dx²) / simplicité / coût O(N)

**Pratique** : Tests validation passent, aucun problème critique identifié

#### Amélioration recommandée

**Ajouter ordre 4 optionnel** (via paramètre `parameters.yaml`) :

```yaml
numerical_parameters:
  derivative_order: 2  # Ou 4 si précision critique
```

**Effort** : ~2-3h implémentation + tests

#### Impact sur autres décisions

**Aucun impact** : D2 déjà résolu, autres décisions compatibles

**Priorité** : 🟢 **BASSE** (amélioration optionnelle)

***

## D3 : Gestion bords grille spatiale

### 1. Contexte physique et numérique

#### Problème posé

Fonction d'onde ψ(x,t) discrétisée sur grille finie `[x_min, x_max]`.

**Question** : Que valent ψ(x_min), ψ(x_max), et leurs dérivées ?

Le cours ne spécifie **aucune condition aux bords** explicitement.

#### Contraintes physiques

1. **Normalisation** : `∫|ψ|² dx = 1`
    - Intégrale sur ℝ entier en théorie
    - Grille finie : approximation
2. **Conservation probabilité** : Flux sortant minimal
    - Si ψ non nul aux bords : probabilité "fuit"
3. **Compatibilité potentiel** :
    - Puits infini : ψ(bords) = 0 forcément
    - Particule libre : ψ(bords) petit si grille assez grande

### 2. État actuel de l'implémentation

**Status** : ✅ **RÉSOLU** (journal, section 2, D3)

#### Implémenté

✅ **Conditions Dirichlet par défaut** :

```
ψ(x_min) = ψ(x_max) = 0
```

✅ **Padding dans fonctions gradient/laplacien** (`utils/numerical.py`)

- Évite accès hors grille (ψᵢ₋₁ quand i=0)

✅ **Documentation ajoutée** dans docstrings

#### Justification (journal)

- Simplicité implémentation
- Compatible puits infini (système déjà implémenté)
- Particule libre : grille suffisamment grande → ψ≈0 aux bords naturellement


### 3. Options techniques possibles

#### Option 1 : Conditions Dirichlet (IMPLÉMENTÉ)

**Description** : `ψ(x_min) = ψ(x_max) = 0`

**Interprétation physique** : Murs impénétrables (puits infini implicite)

**Avantages** :

- ✅ Simplicité maximale
- ✅ Compatible puits infini (Complément HI cours)
- ✅ Évite flux sortant

**Inconvénients** :

- ⚠️ Réflexions non physiques si ψ atteint bords
    - Particule libre : nécessite grille large (±5σ)

**Compatibilité cours** : ✅ **Compatible puits infini**

**Cas d'usage** :

- ✅ Puits infini
- ✅ Particule libre (si grille grande)
- ❌ États étendus atteignant bords

***

#### Option 2 : Conditions périodiques

**Description** : `ψ(x_min) = ψ(x_max)`, `∂ψ/∂x(x_min) = ∂ψ/∂x(x_max)`

**Interprétation physique** : Topologie circulaire (anneau)

**Avantages** :

- Compatible FFT (nécessaire Option 3, D1)
- Pas de réflexions

**Inconvénients** :

- ❌ **Non physique** pour la plupart des systèmes
    - Particule "revient" par bord opposé
- ❌ Incompatible puits infini

**Compatibilité cours** : ❌ Hypothèse non justifiée

**Verdict** : Réservé split-operator (si implémenté)

***

#### Option 3 : Conditions absorbantes (PML)

**Description** : Couches absorbantes (Perfectly Matched Layers)

- Potentiel imaginaire pur aux bords : `V → V - iΓ(x)`
- Absorbe onde sortante sans réflexion

**Avantages** :

- ✅ Élimine réflexions non physiques
- ✅ Physiquement réaliste (grille = fenêtre sur ℝ)

**Inconvénients** :

- ❌ **Complexité implémentation** élevée
- ❌ Hamiltonien non hermitique (V complexe)
    - ⚠️ Viole Règle R4.5 localement
- ❌ Norme décroît (absorption)

**Compatibilité cours** : ❌ Hamiltonien non hermitique non traité

**Verdict** : **Hors périmètre** - Trop avancé

### 4. Évaluation critique

#### Options à exclure

**PML (Option 3)** : ❌ Complexité + hamiltonien non hermitique

#### Options acceptables

1. **Dirichlet (Option 1)** : ✅ **IMPLÉMENTÉ** - Optimal contexte actuel
2. **Périodiques (Option 2)** : ⚠️ Réservé FFT (si nécessaire)

#### Option actuelle cohérente

**Conditions Dirichlet** : ✅ **CONFIRMÉ COMME OPTIMAL**

### 5. Recommandation

#### Décision recommandée

**Conserver Dirichlet** avec **validation grille adaptée**

#### Justification

- Simplicité maximale
- Compatible tous systèmes actuels (libre, puits infini)
- Tests passent (journal)


#### Amélioration recommandée

**Valider dimensionnement grille** :

```python
def validate_grid_size(psi: np.ndarray, threshold=1e-4):
    """Vérifie ψ négligeable aux bords."""
    assert abs(psi[^0]) < threshold, "ψ(x_min) trop grand ! Agrandir grille."
    assert abs(psi[-1]) < threshold, "ψ(x_max) trop grand ! Agrandir grille."
```

**Règle dimensionnement** (journal, section 9.2) :

- Gaussienne σₓ : grille ±5σₓ (couvre 99.9999%)


#### Impact sur autres décisions

**Aucun impact** : D3 résolu, compatible D1 (Crank-Nicolson), D2 (différences finies)

**Priorité** : 🟢 **RÉSOLUE**

***

## D4 : Construction état fondamental oscillateur harmonique

### 1. Contexte physique et numérique

#### Problème physique posé

Oscillateur harmonique 1D (Chapitre V cours) :

**Hamiltonien** (Règle R6.1) :

```
H = P²/2m + (1/2)mω²X²
```

**Spectre** : `Eₙ = ℏω(n + 1/2)`, n = 0, 1, 2, ...

**États propres |n⟩** définis algébriquement (Règles R6.2, R6.3) :

- État fondamental : `a|0⟩ = 0`
- États excités : `|n⟩ = (a†)ⁿ/√(n!) |0⟩`
- Opérateurs échelle : `a|n⟩ = √n|n-1⟩`, `a†|n⟩ = √(n+1)|n+1⟩`

**Source** : [file:1, Chapitre V, § B-C]

**Question** : Comment obtenir ψₙ(x) = ⟨x|n⟩ en représentation position ?

#### Équations concernées

Le cours **mentionne** (Complément BV) polynômes d'Hermite mais **détails absents** des extraits fournis (Limite L2).

**Formule attendue** (non dans extraits) :

```
ψₙ(x) = (mω/πℏ)^(1/4) · 1/√(2ⁿn!) · Hₙ(√(mω/ℏ)x) · exp(-mωx²/2ℏ)
```

où Hₙ = polynômes Hermite

### 2. État actuel de l'implémentation

**Status** : ✅ **RÉSOLU (choix alternatif)** (journal, section 1.2)

#### Implémenté

✅ **Algèbre opérateurs échelle en base abstraite** :

- `energy_eigenvalue(n)` : Eₙ = ℏω(n+1/2)
- Opérateurs a, a† : Relations commutation, action sur |n⟩
- **Base de Fock** {|n⟩} : Représentation matricielle

✅ **Décision adoptée** : **Option 3 - Base abstraite** (journal)

#### Limites

❌ **Fonctions d'onde ψₙ(x) NON implémentées**

- Polynômes Hermite absents
- **Impact** :
    - ✅ Spectroscopie HO fonctionnelle (niveaux, transitions)
    - ❌ Visualisation ψₙ(x) impossible
    - ❌ Évolution paquets HO en représentation position bloquée


### 3. Options techniques possibles

#### Option 1 : Résolution numérique équation Schrödinger stationnaire

**Description** : Diagonaliser H en représentation position

**Méthode** :

1. Construire matrice H sur grille x
2. Diagonaliser : `H|ψₙ⟩ = Eₙ|ψₙ⟩`
3. Vecteurs propres = ψₙ(x) discrétisés

**Avantages** :

- ✅ Pas besoin formules analytiques
- ✅ Généralisable autres potentiels

**Inconvénients** :

- Coût diagonalisation : O(N³) si matrice pleine, O(N²) si sparse
- Erreurs discrétisation sur ψₙ
- **Choix grille critique** : doit couvrir ψₙ (n grand → étendu)

**Compatibilité cours** : ✅ Résout équation Schrödinger stationnaire (cours Chapitre III, §D-1)

**Coût** : O(N²-N³) selon méthode

***

#### Option 2 : Utiliser `scipy.special.hermite` (externe au cours)

**Description** : Implémenter formule explicite avec bibliothèque Python

```python
from scipy.special import hermite
from numpy.polynomial.hermite import hermval

def psi_n(x, n, m, omega, hbar):
    alpha = np.sqrt(m * omega / hbar)
    norm = (m*omega/(np.pi*hbar))**0.25 / np.sqrt(2**n * factorial(n))
    Hn = hermval(alpha * x, [^0]*n + [^1])  # Polynôme Hermite
    return norm * Hn * np.exp(-alpha**2 * x**2 / 2)
```

**Avantages** :

- ✅ Précision analytique (machine)
- ✅ Efficace
- ✅ Lève Limite L2

**Inconvénients** :

- ❌ **Source externe au cours**
    - Polynômes Hermite mentionnés mais non détaillés
    - Bibliothèque scipy = connaissance externe

**Compatibilité cours** : ⚠️ **Formules non fournies** dans extraits

**Verdict** : Acceptable si **explicitement documenté** comme extension

***

#### Option 3 : Base abstraite {|n⟩} (IMPLÉMENTÉ)

**Description** : Travailler uniquement avec matrices

**Représentation** :

- |0⟩ = [1, 0, 0, ..., 0]ᵀ
- |1⟩ = [0, 1, 0, ..., 0]ᵀ
- ...
- Troncature `n_max_fock` (journal : 50)

**Opérateurs** :

```
a = matrice tridiagonale inférieure (éléments √n)
a† = matrice tridiagonale supérieure (éléments √(n+1))
H = ℏω(a†a + 1/2·I)
```

**Avantages** :

- ✅ **Strictement basé sur cours** (Chapitre V)
- ✅ Algèbre exacte (Règles R6.2, R6.3)
- ✅ Suffit pour spectroscopie

**Inconvénients** :

- ❌ Pas de ψₙ(x) → pas de visualisation
- ❌ Évolution temporelle limitée (états stationnaires seulement)

**Compatibilité cours** : ✅ **PARFAITE** - Utilise uniquement définitions algébriques

***

### 4. Évaluation critique

#### Options à exclure

**Aucune** - Toutes valides selon contraintes

#### Options acceptables

1. **Base abstraite (Option 3)** : ✅ **IMPLÉMENTÉ** - Cohérent avec cours fourni
2. **Scipy Hermite (Option 2)** : ✅ **Extension future** - Si visualisation nécessaire
3. **Diagonalisation (Option 1)** : ⚠️ Acceptable mais coûteux

#### Option actuelle cohérente

**Base abstraite** : ✅ **CONFIRMÉ COMME OPTIMAL** dans contexte actuel

**Raisons** :

- Respecte strictement extraits cours fournis (pas d'ajout externe)
- Algèbre opérateurs échelle = cœur physique HO
- Suffisant pour applications actuelles (spectres, transitions)


### 5. Recommandation

#### Décision recommandée

**Conserver base abstraite** + **Extension Option 2 future** si visualisation requise

#### Justification synthétique

**Physique** : Algèbre échelle = contenu essentiel Chapitre V (complet dans implémentation)

**Cohérence cours** : Pas d'ajout non justifié (Hermite absents extraits)

**Pratique** : Fonctionnel pour applications actuelles

#### Extension future recommandée (journal, section 5.1, priorité 1)

**Si visualisation ψₙ(x) nécessaire** :

```python
# Ajouter dans systems/harmonic_oscillator.py
def wavefunction_position(self, n: int, x_grid: np.ndarray) -> np.ndarray:
    """
    Fonction d'onde ψₙ(x) en représentation position.
    
    EXTENSION : Utilise scipy.special.hermite (hors extraits cours).
    Formule complète dans Complément BV (Tome I).
    """
    from scipy.special import eval_hermite
    # ... implémentation
```

**Documentation obligatoire** : Préciser source externe

**Effort** : ~1-2h implémentation + tests

#### Impact sur autres décisions

**D1 (évolution)** : ⚠️ Attention

- Actuellement : évolution via phases `exp(-iEₙt/ℏ)` (états stationnaires)
- Avec ψₙ(x) : pourrait évoluer paquets HO en représentation position
    - Nécessite Crank-Nicolson opérationnel

**Autres** : Aucun impact

**Priorité** : 🟡 **MOYENNE** (amélioration optionnelle)

***

## D5 : Tirage aléatoire mesures

### 1. Contexte physique et numérique

#### Problème posé

Lors d'une **mesure quantique** (4ème postulat), le résultat est aléatoire.

**Probabilités** (Règle R2.2) :

```
P(aₙ) = |⟨uₙ|ψ⟩|²
```

où aₙ = valeurs propres observable A, uₙ = vecteurs propres

**Source** : [file:1, Chapitre III, § B-3-b]

**Question implémentation** : Comment simuler ce tirage aléatoire ?

#### Contraintes physiques

1. **Respect probabilités quantiques** : Tirage selon distribution P(aₙ)
2. **Normalisation** : `Σ P(aₙ) = 1` (validation pré-tirage)
3. **Réduction paquet d'ondes** (Règle R2.3) : Après mesure aₙ, état devient projection

**Le cours ne spécifie PAS** méthode de tirage (aspect implémentation pure).

### 2. État actuel de l'implémentation

**Status** : ✅ **RÉSOLU** (journal, section 1.1, D5)

#### Implémenté

✅ **`np.random.choice` avec seed optionnel** :

```python
def measure_once(self, state: QuantumState, random_seed=None) -> tuple:
    """
    Simule UNE mesure selon probabilités quantiques.
    Retourne (valeur_mesurée, état_après_mesure).
    """
    if random_seed is not None:
        np.random.seed(random_seed)
    
    probabilities = self.compute_probabilities(state)
    eigenvalues = list(probabilities.keys())
    probs = list(probabilities.values())
    
    measured_value = np.random.choice(eigenvalues, p=probs)
    state_after = self.apply_reduction(state, measured_value)
    
    return measured_value, state_after
```

✅ **Logging complet** (journal, amélioration) :

- Séquence mesures enregistrée
- Statistiques finales (fréquences observées vs théoriques)

✅ **Tests validation** (journal, section 1.1) :

- Test χ² : distribution empirique vs théorique
- Réduction paquet : mesures successives donnent même résultat


#### Justification (journal)

- Seed optionnel → reproductibilité (tests, debug)
- `np.random.choice` : implémentation standard Python
- Validation statistique : test χ² sur N=1000 mesures


### 3. Options techniques possibles

#### Option 1 : `np.random.choice` (IMPLÉMENTÉ)

**Description** : Fonction standard NumPy pour tirage discret

```python
np.random.choice(values, p=probabilities)
```

**Avantages** :

- ✅ Simplicité maximale
- ✅ Bibliothèque standard (NumPy)
- ✅ Optimisée (C interne)
- ✅ Seed pour reproductibilité

**Inconvénients** :

- Aucun (pour cas d'usage)

**Compatibilité cours** : ✅ Implémente postulat mesure (Règle R2.2)

**Verdict** : ✅ **OPTIMAL**

***

#### Option 2 : Méthode inverse (manuelle)

**Description** : Implémenter tirage manuellement

```python
def sample(eigenvalues, probabilities):
    cumulative = np.cumsum(probabilities)
    u = np.random.uniform(0, 1)
    index = np.searchsorted(cumulative, u)
    return eigenvalues[index]
```

**Avantages** :

- Contrôle total algorithme
- Pédagogique

**Inconvénients** :

- Réinvente roue (np.random.choice fait déjà ça)
- Pas de gain

**Verdict** : Inutile

***

#### Option 3 : Générateur quantique "vrai" (hardware)

**Description** : Utiliser générateur aléatoire quantique matériel

**Problème** : Hors périmètre simulation logicielle

**Verdict** : Non applicable

### 4. Évaluation critique

#### Option la plus cohérente

**`np.random.choice`** : ✅ **IMPLÉMENTÉ ET OPTIMAL**

**Raisons** :

- Implémentation standard robuste
- Tests validation passent (χ², réduction)
- Seed optionnel répond besoin reproductibilité


### 5. Recommandation

#### Décision recommandée

**Conserver `np.random.choice`** - Aucune modification nécessaire

#### Justification

- Implémentation correcte physiquement (Règle R2.2)
- Tests validation réussis (journal, section 6.1)
- Simplicité + robustesse


#### Point d'attention : Seed en production

**Question ouverte** : Faut-il exposer `random_seed` dans `parameters.yaml` ?

**Recommandation** : ❌ **NON**

- Seed = paramètre debug/tests uniquement
- Mesures quantiques = intrinsèquement aléatoires
- En production : laisser aléatoire (seed=None par défaut)

**Si besoin reproductibilité** : Passer seed explicitement dans code expérience

```python
# Dans script expérience
measurement = QuantumMeasurement(observable)
result, state_after = measurement.measure_once(state, random_seed=42)  # Debug
```


#### Impact sur autres décisions

**Aucun** - D5 résolu, indépendant des autres

**Priorité** : 🟢 **RÉSOLUE**

***

## Tableau de synthèse final

| Décision | Problème | Option retenue | Statut | Priorité | Effort estimé | Impact autres décisions |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| **D1** | Schéma intégration temporelle | **Crank-Nicolson** | ⚠️ **À IMPLÉMENTER** | 🔴 **HAUTE** | 6-8h | Compatible D2, D3 ; améliore validation |
| **D2** | Gradient/laplacien | Différences finies ordre 2 (+ ordre 4 optionnel) | ✅ **RÉSOLU** | 🟢 Basse (amélioration) | 2-3h (ordre 4) | Aucun |
| **D3** | Bords grille spatiale | Conditions Dirichlet | ✅ **RÉSOLU** | 🟢 Résolue | - | Compatible D1, D2 |
| **D4** | État fondamental HO | Base abstraite { | n⟩} (+ scipy.hermite futur) | ✅ **RÉSOLU** | 🟡 Moyenne (extension) | 1-2h (ψₙ(x)) |
| **D5** | Tirage aléatoire mesures | `np.random.choice` | ✅ **RÉSOLU** | 🟢 Résolue | - | Aucun |

### Dépendances entre décisions

```
D1 (Crank-Nicolson) ─── bloque ──→ Évolution continue générale
  │
  ├─ nécessite → D2 (laplacien)  ✅ Résolu
  ├─ nécessite → D3 (bords)      ✅ Résolu
  └─ améliore  → Validation conservation (95%→100%)

D4 (HO) ─── si ψₙ(x) implémenté ──→ nécessite D1 pour évolution paquets

D2, D3, D5 : Indépendants et résolus
```


### Actions immédiates recommandées

**Priorité 1** : 🔴 **Implémenter Crank-Nicolson (D1)**

- Bloqueur critique identifié (journal, section 1.2)
- Toutes dépendances résolues (D2, D3)
- Effort : ~1 semaine développement + validation

**Priorité 2** : 🟡 **Polynômes Hermite (D4 extension)**

- Après D1 opérationnel
- Permet visualisations ψₙ(x)
- Effort : ~2-3h

**Priorité 3** : 🟢 **Différences finies ordre 4 (D2 amélioration)**

- Optionnel, si précision critique
- Effort : ~2-3h

***

## Conclusion

### Décisions validées

✅ **D2, D3, D5** : Implémentations actuelles **CONFIRMÉES COMME OPTIMALES**

- Testées et fonctionnelles (journal, section 6)
- Compatibles cours et architecture
- Aucune modification nécessaire


### Décision critique

⚠️ **D1** : **Action requise urgente**

- Crank-Nicolson = **seule option physiquement acceptable**
    - Conservation norme exacte (Règle R5.1)
    - Stabilité inconditionnelle
    - Compatible architecture actuelle
- Bloque évolution continue (fonctionnalité centrale)


### Décision optionnelle

🟡 **D4 extension** : Polynômes Hermite

- Base abstraite actuelle = cohérente avec cours fourni
- Extension justifiable si visualisation nécessaire
- Dépend de D1 pour pleine utilité

***

**Document produit le** : 17 décembre 2025
**Basé exclusivement sur** : Cours Cohen-Tannoudji + Journal implémentation
**Validation** : Toutes options traçables aux contraintes physiques du cours
<span style="display:none">[^2]</span>

<div align="center">⁂</div>

[^1]: Journal-des-changements-et-ameliorations.md

[^2]: Mecanique-quantique.-Tome-I.pdf

