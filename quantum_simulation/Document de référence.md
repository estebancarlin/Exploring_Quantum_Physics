# Document de référence unique : Implémentation de simulations de mécanique quantique

## Basé sur Cohen-Tannoudji, Diu, Laloë - Mécanique Quantique Tome I

**Version** : 1.0
**Date** : 16 décembre 2025
**Statut** : Document de référence pour implémentation Python

***

## 1. Cadre théorique issu du cours

### 1.1 Principes fondamentaux

#### 1.1.1 États quantiques

**Source** : [file:1, Chapitre III, § B-1]

Un système quantique est décrit par un **vecteur d'état** noté `|ψ⟩` (notation de Dirac), appartenant à un espace de Hilbert E appelé **espace des états**.

**Propriétés essentielles** :

- Linéarité : si `|ψ₁⟩` et `|ψ₂⟩` sont états possibles, alors `α|ψ₁⟩ + β|ψ₂⟩` aussi (α, β complexes)
- Normalisation : `⟨ψ|ψ⟩ = 1` pour un état physique (convention)
- **Équivalence physique** : `|ψ⟩` et `e^(iθ)|ψ⟩` décrivent le même état (phase globale sans effet) - Source : [file:1, Chapitre III, § B-3-b-γ]

**Représentations possibles** :

- **Représentation position** {|r⟩} : l'état s'écrit comme fonction d'onde `ψ(r) = ⟨r|ψ⟩` - Source : [file:1, Chapitre II, § E-1]
- **Représentation impulsion** {|p⟩} : `φ(p) = ⟨p|ψ⟩` - Source : [file:1, Chapitre II, § E-2]
- **Base d'états propres** : décomposition `|ψ⟩ = Σ cₙ|uₙ⟩` sur états propres d'une observable

**Produit scalaire** : `⟨φ|ψ⟩` (complexe) - Source : [file:1, Chapitre II, § B-2-c]

- Propriété : `⟨φ|ψ⟩* = ⟨ψ|φ⟩`
- Norme : `||ψ|| = √⟨ψ|ψ⟩`


#### 1.1.2 Observables

**Source** : [file:1, Chapitre II, § D-1]

Une **observable** A est un opérateur hermitique (A† = A) agissant sur l'espace des états.

**Propriétés garanties** :

- Valeurs propres réelles : `A|uₙ⟩ = aₙ|uₙ⟩` avec aₙ ∈ ℝ
- Vecteurs propres orthonormés : `⟨uᵢ|uⱼ⟩ = δᵢⱼ`
- Les {|uₙ⟩} forment une base complète de E (relation de fermeture : `Σ|uₙ⟩⟨uₙ| = 𝟙`)

**Observables fondamentales** :

1. **Position** R = (X, Y, Z) - Source : [file:1, Chapitre II, § E-1]
2. **Impulsion** P = (Pₓ, Pᵧ, Pᵤ) - Source : [file:1, Chapitre II, § E-2]
    - En représentation position : `P = -iℏ∇`
3. **Hamiltonien** H (énergie totale) - Source : [file:1, Chapitre III, § B-4]

**Relations de commutation canoniques** - Source : [file:1, Chapitre III, § B-5-a] :

```
[Rᵢ, Rⱼ] = 0
[Pᵢ, Pⱼ] = 0
[Rᵢ, Pⱼ] = iℏδᵢⱼ
```

**Compatibilité** : A et B compatibles ⟺ `[A,B] = 0` ⟺ bases propres communes possibles - Source : [file:1, Chapitre III, § C-6-a]

#### 1.1.3 Mesure quantique

**4ème Postulat (probabilités)** - Source : [file:1, Chapitre III, § B-3-b]

Pour un système dans l'état normé `|ψ⟩`, la mesure de l'observable A donne :

**Cas spectre discret non dégénéré** :

- Valeurs possibles : valeurs propres {aₙ} de A
- Probabilité d'obtenir aₙ : `P(aₙ) = |⟨uₙ|ψ⟩|²` où A|uₙ⟩ = aₙ|uₙ⟩

**Cas spectre continu** :

- Densité de probabilité : `dP/dα = |⟨vₐ|ψ⟩|²`

**5ème Postulat (réduction du paquet d'ondes)** - Source : [file:1, Chapitre III, § B-3-c]

Après mesure donnant aₙ, l'état devient immédiatement :

```
|ψ'⟩ = Pₙ|ψ⟩ / √⟨ψ|Pₙ|ψ⟩
```

où Pₙ est le projecteur sur le sous-espace propre associé à aₙ.

**Conséquence** : seconde mesure immédiate de A redonne aₙ avec certitude.

#### 1.1.4 Évolution temporelle

**6ème Postulat (équation de Schrödinger)** - Source : [file:1, Chapitre III, § B-4]

L'évolution de l'état est régie par :

```
iℏ d|ψ(t)⟩/dt = H(t)|ψ(t)⟩
```

**Cas conservatif** (H indépendant de t) - Source : [file:1, Chapitre III, § D-2]

Méthode de résolution par décomposition spectrale :

1. Diagonaliser H : trouver {Eₙ, |φₙ⟩} tels que H|φₙ⟩ = Eₙ|φₙ⟩
2. Décomposer état initial : `|ψ(t₀)⟩ = Σ cₙ(t₀)|φₙ⟩`
3. Évolution : `cₙ(t) = cₙ(t₀) exp(-iEₙ(t-t₀)/ℏ)`
4. État au temps t : `|ψ(t)⟩ = Σ cₙ(t)|φₙ⟩`

**États stationnaires** : états propres de H ne changent que par phase globale (physiquement invariants) - Source : [file:1, Chapitre III, § D-2-b]

#### 1.1.5 Grandeurs moyennes et incertitudes

**Valeur moyenne** - Source : [file:1, Chapitre III, § C-4]

```
⟨A⟩ = ⟨ψ|A|ψ⟩
```

**Écart quadratique moyen (incertitude)** - Source : [file:1, Chapitre III, § C-5]

```
ΔA = √(⟨A²⟩ - ⟨A⟩²)
```

**Relations d'incertitude de Heisenberg** - Source : [file:1, Chapitre III, § C-5]

```
ΔX · ΔPₓ ≥ ℏ/2
ΔY · ΔPᵧ ≥ ℏ/2
ΔZ · ΔPᵤ ≥ ℏ/2
```

**Théorème d'Ehrenfest** - Source : [file:1, Chapitre III, § D-1-d]

```
d⟨R⟩/dt = ⟨P⟩/m
d⟨P⟩/dt = -⟨∇V(R)⟩
```

Lien avec mécanique classique : valeurs moyennes obéissent à équations classiques.

#### 1.1.6 Conservation de la probabilité

**Équation de continuité** - Source : [file:1, Chapitre III, § D-1-c]

```
∂ρ/∂t + ∇·J = 0
```

où :

- `ρ(r,t) = |ψ(r,t)|²` (densité de probabilité)
- `J(r,t) = (ℏ/2mi)[ψ*∇ψ - ψ∇ψ*] = (1/m)Re(ψ* (ℏ/i)∇ψ)` (courant de probabilité)

**Conséquence** : `∫ρ(r,t) d³r = constante = 1` (norme conservée)

### 1.2 Systèmes physiques spécifiques couverts

#### 1.2.1 Particule libre (V = 0)

**Hamiltonien** - Source : [file:1, Chapitre I, § C]

```
H = P²/2m
```

En représentation position :

```
H = -ℏ²/2m Δ
```

**États propres** : ondes planes `exp(ik·r)` de valeur propre `E = ℏ²k²/2m`

#### 1.2.2 Particule dans potentiel scalaire V(r)

**Hamiltonien** - Source : [file:1, Chapitre III, § B-5-b]

```
H = P²/2m + V(R)
```

En représentation position :

```
Hψ = -ℏ²/2m Δψ + V(r)ψ
```

**Équation de Schrödinger (représentation position)** - Source : [file:1, Chapitre I, § B-2]

```
iℏ ∂ψ(r,t)/∂t = -ℏ²/2m Δψ(r,t) + V(r,t)ψ(r,t)
```


#### 1.2.3 Oscillateur harmonique à 1D

**Source** : [file:1, Chapitre V]

**Hamiltonien** :

```
H = P²/2m + (1/2)mω²X²
```

**Opérateurs d'échelle (création/annihilation)** :

```
a = √(mω/2ℏ)(X + i/(mω)P)
a† = √(mω/2ℏ)(X - i/(mω)P)
```

**Relation de commutation** : `[a, a†] = 1`

**Réécriture hamiltonien** : `H = ℏω(N + 1/2)` où `N = a†a`

**Spectre d'énergie** :

```
Eₙ = ℏω(n + 1/2),  n = 0, 1, 2, ...
```

**États propres |n⟩** :

- État fondamental : `a|0⟩ = 0`
- États excités : `|n⟩ = (a†)ⁿ/√(n!) |0⟩`
- Action des opérateurs échelle :

```
a|n⟩ = √n |n-1⟩
a†|n⟩ = √(n+1) |n+1⟩
```


#### 1.2.4 Atome d'hydrogène

**Source** : [file:1, Chapitre VII] - **PARTIEL dans extraits fournis**

**Résultats disponibles** (Complément CI) :

- Rayon de Bohr : `a₀ = ℏ²/(me²)` où `e² = q²/(4πε₀)`
- Énergie fondamentale : `E₀ = -me⁴/(2ℏ²)`

**LIMITE** : Méthode de résolution complète, fonctions d'onde explicites, spectre complet non fournis dans extraits.

### 1.3 Notations et conventions retenues

#### Constantes physiques

- h : constante de Planck (J·s) - ordre grandeur : 6.62×10⁻³⁴
- ℏ = h/(2π) : constante réduite
- m : masse particule (kg)
- ω : pulsation (rad/s)
- k : vecteur d'onde (m⁻¹)
- E : énergie (J)
- p : impulsion (kg·m/s)


#### Opérateurs et états

- `|ψ⟩` : ket (vecteur d'état)
- `⟨ψ|` : bra (dual)
- `⟨φ|ψ⟩` : produit scalaire
- `A, B, H` : opérateurs (majuscules)
- `[A,B] = AB - BA` : commutateur
- `𝟙` : opérateur identité
- `δᵢⱼ` : symbole de Kronecker
- `δ(x)` : distribution de Dirac


#### Représentations

- `ψ(r,t)` : fonction d'onde en représentation position
- `φ(p,t)` : fonction d'onde en représentation impulsion
- `r = (x, y, z)` : vecteur position
- `∇ = (∂/∂x, ∂/∂y, ∂/∂z)` : gradient
- `Δ = ∇² = ∂²/∂x² + ∂²/∂y² + ∂²/∂z²` : laplacien

***

## 2. Règles physiques implémentables

### 2.1 Relations fondamentales

#### Règle R1.1 : Relations de Planck-Einstein

**Énoncé** : Correspondance onde-corpuscule pour photons
**Formulation** :

```
E = ℏω
p = ℏk  (vectoriel)
```

**Source** : [file:1, Chapitre I, § A-1]
**Contrainte numérique** : Cohérence dimensionnelle h en J·s

#### Règle R1.2 : Relations de Louis de Broglie

**Énoncé** : Onde de matière associée à particule
**Formulation** :

```
λ = h/p
ν = E/h
```

**Source** : [file:1, Chapitre I, § B-1]
**Usage** : Conversions entre paramètres ondulatoires/corpusculaires

#### Règle R1.3 : Relations de commutation canoniques

**Énoncé** : Structure algébrique position-impulsion
**Formulation** :

```
[Rᵢ, Rⱼ] = 0
[Pᵢ, Pⱼ] = 0
[Rᵢ, Pⱼ] = iℏδᵢⱼ
```

**Source** : [file:1, Chapitre III, § B-5-a]
**Contrainte** : Vérifier numériquement `||[Rᵢ,Pⱼ] - iℏδᵢⱼ|| < ε` sur états tests
**Invariant** : Structure préservée par toute transformation unitaire

### 2.2 Interprétation probabiliste

#### Règle R2.1 : Densité de probabilité

**Énoncé** : Probabilité de présence en représentation position
**Formulation** :

```
ρ(r,t) = |ψ(r,t)|²
P(r ∈ V, t) = ∫_V ρ(r,t) d³r
```

**Source** : [file:1, Chapitre I, § B-2]
**Contrainte** : Normalisation `∫ρ d³r = 1`

#### Règle R2.2 : Probabilité de mesure (spectre discret)

**Énoncé** : Résultat mesure observable à spectre discret
**Formulation** :

```
P(aₙ) = |⟨uₙ|ψ⟩|²    (cas non dégénéré)
P(aₙ) = Σᵢ |⟨uᵢ⁽ⁿ⁾|ψ⟩|²  (cas dégénéré, somme sur sous-espace propre)
```

**Source** : [file:1, Chapitre III, § B-3-b]
**Contrainte** : `Σₙ P(aₙ) = 1` (vérification post-calcul)

#### Règle R2.3 : Réduction du paquet d'ondes

**Énoncé** : Modification état après mesure
**Formulation** :

```
|ψ⟩  --mesure donne aₙ-->  |ψ'⟩ = Pₙ|ψ⟩ / √⟨ψ|Pₙ|ψ⟩
```

où Pₙ = projecteur sur sous-espace propre(aₙ)
**Source** : [file:1, Chapitre III, § B-3-c]
**Contrainte** : `⟨ψ'|ψ'⟩ = 1` (normalisation)
**Effet** : Mesure répétée immédiate donne aₙ avec probabilité 1

### 2.3 Évolution temporelle

#### Règle R3.1 : Équation de Schrödinger (forme abstraite)

**Énoncé** : Évolution unitaire
**Formulation** :

```
iℏ d|ψ(t)⟩/dt = H(t)|ψ(t)⟩
```

**Source** : [file:1, Chapitre III, § B-4]
**Contrainte** : Conservation norme `d⟨ψ|ψ⟩/dt = 0` (découle hermiticité H)

#### Règle R3.2 : Équation de Schrödinger (représentation position)

**Énoncé** : Forme différentielle
**Formulation** :

```
iℏ ∂ψ(r,t)/∂t = [-ℏ²/2m Δ + V(r,t)] ψ(r,t)
```

**Source** : [file:1, Chapitre I, § B-2]
**Contrainte numérique** :

- Discrétisation Δ (différences finies, FFT, etc.) - **choix implémentation**
- Schéma intégration temporelle (Euler, RK, split-operator) - **choix implémentation**


#### Règle R3.3 : Évolution par décomposition spectrale

**Énoncé** : Solution formelle système conservatif
**Formulation** :

```
|ψ(t₀)⟩ = Σ cₙ(t₀)|φₙ⟩
⇒ |ψ(t)⟩ = Σ cₙ(t₀) exp[-iEₙ(t-t₀)/ℏ] |φₙ⟩
```

où H|φₙ⟩ = Eₙ|φₙ⟩
**Source** : [file:1, Chapitre III, § D-2-a]
**Contrainte** : Nécessite diagonalisation préalable de H

#### Règle R3.4 : États stationnaires

**Énoncé** : États propres H physiquement invariants
**Formulation** :

```
|ψ(t)⟩ = exp[-iE(t-t₀)/ℏ] |φ⟩  où H|φ⟩ = E|φ⟩
```

Phase globale → toutes observables constantes
**Source** : [file:1, Chapitre III, § D-2-b]
**Test** : `⟨A⟩(t) = ⟨A⟩(t₀)` pour toute observable A

### 2.4 Observables et mesure

#### Règle R4.1 : Valeur moyenne

**Énoncé** : Espérance mathématique observable
**Formulation** :

```
⟨A⟩ = ⟨ψ|A|ψ⟩
```

État normé requis.
**Source** : [file:1, Chapitre III, § C-4]

#### Règle R4.2 : Écart quadratique moyen

**Énoncé** : Incertitude quantique
**Formulation** :

```
ΔA = √(⟨A²⟩ - ⟨A⟩²)
```

**Source** : [file:1, Chapitre III, § C-5]
**Contrainte** : ΔA ≥ 0 par construction

#### Règle R4.3 : Relations d'incertitude de Heisenberg

**Énoncé** : Limites fondamentales mesure simultanée
**Formulation** :

```
ΔX · ΔPₓ ≥ ℏ/2
ΔY · ΔPᵧ ≥ ℏ/2
ΔZ · ΔPᵤ ≥ ℏ/2
```

**Source** : [file:1, Chapitre III, § C-5]
**Test validation** : Calculer ΔX, ΔPₓ sur tout état, vérifier inégalité

#### Règle R4.4 : Théorème d'Ehrenfest

**Énoncé** : Équations classiques pour valeurs moyennes
**Formulation** :

```
d⟨R⟩/dt = ⟨P⟩/m
d⟨P⟩/dt = -⟨∇V(R)⟩
```

**Source** : [file:1, Chapitre III, § D-1-d]
**Test validation** : Calculer dérivées temporelles numériquement, comparer membres

#### Règle R4.5 : Hermiticité observables

**Énoncé** : Contrainte structure opérateurs mesurables
**Formulation** :

```
A† = A  ⟺  ⟨φ|A|ψ⟩* = ⟨ψ|A|φ⟩
```

**Source** : [file:1, Chapitre II, § D-1]
**Test** : Vérifier `||A† - A|| < ε` sur base discrète

### 2.5 Conservation

#### Règle R5.1 : Conservation de la probabilité

**Énoncé** : Norme état constante
**Formulation** :

```
d⟨ψ(t)|ψ(t)⟩/dt = 0
```

**Source** : [file:1, Chapitre III, § D-1-c]
**Conséquence** : `∫|ψ(r,t)|² d³r = constante`

#### Règle R5.2 : Équation de continuité

**Énoncé** : Conservation locale probabilité
**Formulation** :

```
∂ρ/∂t + ∇·J = 0
```

avec :

```
ρ(r,t) = |ψ(r,t)|²
J(r,t) = (ℏ/2mi)[ψ*∇ψ - ψ∇ψ*]
```

**Source** : [file:1, Chapitre III, § D-1-c]
**Test numérique** : Calculer membres, vérifier somme ≈ 0

### 2.6 Oscillateur harmonique (règles spécifiques)

#### Règle R6.1 : Spectre oscillateur

**Énoncé** : Quantification énergie
**Formulation** :

```
Eₙ = ℏω(n + 1/2),  n ∈ ℕ
```

**Source** : [file:1, Chapitre V, § B]

#### Règle R6.2 : Algèbre opérateurs échelle

**Énoncé** : Relations définissant a, a†
**Formulation** :

```
[a, a†] = 1
H = ℏω(a†a + 1/2)
```

**Source** : [file:1, Chapitre V, § B]
**Test** : Vérifier commutateur sur base tronquée

#### Règle R6.3 : Action échelle sur états propres

**Énoncé** : Construction récursive états
**Formulation** :

```
a|n⟩ = √n |n-1⟩
a†|n⟩ = √(n+1) |n+1⟩
a|0⟩ = 0
```

**Source** : [file:1, Chapitre V, § C]
**Usage** : Construction |n⟩ = (a†)ⁿ/√(n!) |0⟩

***

## 3. Traduction logicielle des concepts physiques

### 3.1 Correspondances fondamentales

| Concept physique | Module | Classe | Responsabilité |
| :-- | :-- | :-- | :-- |
| État quantique | ψ⟩ | `core.state` | `QuantumState` (abstraite) |
| Fonction d'onde ψ(r) | `core.state` | `WaveFunctionState` | Représentation position sur grille |
| Décomposition Σcₙ | uₙ⟩ | `core.state` | `EigenStateBasis` |
| Observable A | `core.operators` | `Observable` (abstraite) | Application, valeurs propres, ⟨A⟩ |
| Position R | `core.operators` | `PositionOperator` | Multiplication par r |
| Impulsion P | `core.operators` | `MomentumOperator` | -iℏ∇ en représentation position |
| Hamiltonien H | `core.operators` | `Hamiltonian` | P²/2m + V(R), évolution |
| Évolution Schrödinger | `dynamics.evolution` | `TimeEvolution` | Intégration iℏ∂ψ/∂t = Hψ |
| Mesure + réduction | `dynamics.measurement` | `QuantumMeasurement` | Probabilités + projection |
| Oscillateur harmonique | `systems.harmonic_oscillator` | `HarmonicOscillator` | Spectre + opérateurs a, a† |
| Expérience complète | `experiments.base` | `Experiment` | Préparation→évolution→mesure |
| Validation Heisenberg | `validation.heisenberg` | `HeisenbergValidator` | Test ΔX·ΔP ≥ ℏ/2 |
| Conservation | `validation.conservation` | `ConservationValidator` | Test ∂ρ/∂t + ∇·J = 0 |

### 3.2 Hypothèses numériques fondamentales

#### H1 : Discrétisation spatiale

**Décision** : Grille uniforme en représentation position
**Paramètres** : `nx, ny, nz, xmin, xmax, ...` (dans `parameters.yaml`)
**Justification** : Cours donne équations continues, discrétisation = choix implémentation
**Impact** :

- Fonction d'onde → tableau numpy
- Intégrales → sommes discrètes avec poids `dx·dy·dz`
- Dérivées → différences finies (ordre à choisir)


#### H2 : Troncature base de Fock (oscillateur)

**Décision** : {|n⟩, n=0..nₘₐₓ} au lieu de ℕ complet
**Paramètre** : `n_max_fock` (yaml)
**Justification** : États |n> grands : contributions négligeables pour états basse énergie
**Test** : Vérifier contributions |n>nₘₐₓ < seuil tolérance

#### H3 : Schéma d'intégration temporelle

**LIMITE** : Non spécifié dans cours
**Choix possibles** : Euler, Runge-Kutta, split-operator, Crank-Nicolson
**Décision actuelle** : **À implémenter selon stabilité/précision requise**
**Paramètre** : `dt` (yaml), méthode en dur dans code

#### H4 : Représentation opérateurs

**Décision** : Matrices pour espaces dimension finie, fonctions pour espaces continus
**Exemples** :

- Position sur grille : multiplication élément par élément
- Impulsion sur grille : FFT ou différences finies
- Hamiltonien oscillateur : matrice (nₘₐₓ+1)×(nₘₐₓ+1)


#### H5 : Tolérance numérique

**Décision** : Tolérance par défaut = 1e-10 (modifiable yaml)
**Usage** :

- Tests normalisation `|⟨ψ|ψ⟩ - 1| < tol`
- Tests hermiticité `||A† - A|| < tol`
- Validation Heisenberg `ΔX·ΔP - ℏ/2 > -tol`


### 3.3 Séparation responsabilités

#### Couche `core/` : Objets quantiques purs

**Rôle** : Définir structures mathématiques sans méthodes numériques spécifiques
**Ne contient PAS** : Schémas intégration, choix grille, visualisation
**Contient** : Classes abstraites, relations algébriques, invariants physiques
**Exemples** : `QuantumState.inner_product()`, `Observable.is_hermitian()`

#### Couche `dynamics/` : Processus physiques

**Rôle** : Implémenter évolution + mesure selon postulats
**Dépend de** : `core/` (états, opérateurs)
**Ne contient PAS** : Détails systèmes spécifiques (potentiels, etc.)
**Exemples** : `TimeEvolution.evolve_wavefunction()`, `QuantumMeasurement.measure_once()`

#### Couche `systems/` : Systèmes physiques particuliers

**Rôle** : Définir hamiltoniens, états propres pour systèmes du cours
**Dépend de** : `core/`, `dynamics/`
**Exemples** : `HarmonicOscillator`, `PotentialWell`, `FreeParticle`

#### Couche `experiments/` : Orchestration simulations

**Rôle** : Séquencer préparation→évolution→mesure→analyse
**Dépend de** : toutes couches précédentes
**Exemples** : `WavePacketEvolution`, `MeasurementStatistics`

#### Couche `validation/` : Tests physiques

**Rôle** : Vérifier invariants, relations d'incertitude, conservation
**Dépend de** : `core/`, `dynamics/`
**Indépendante de** : systèmes/expériences spécifiques
**Exemples** : `HeisenbergValidator`, `ConservationValidator`

***

## 4. Architecture logicielle globale

### 4.1 Organisation dossiers

```
quantum_simulation/
│
├── config/
│   └── parameters.yaml              # Configuration unique
│
├── core/                             # Fondations quantiques
│   ├── __init__.py
│   ├── constants.py                  # h, ℏ, m, ...
│   ├── state.py                      # QuantumState, WaveFunctionState, EigenStateBasis
│   ├── operators.py                  # Observable, Position, Momentum, Hamiltonian
│   └── hilbert_space.py              # Bases, projections, produits tensoriels
│
├── dynamics/                         # Évolution + mesure
│   ├── __init__.py
│   ├── evolution.py                  # TimeEvolution
│   └── measurement.py                # QuantumMeasurement
│
├── systems/                          # Systèmes physiques
│   ├── __init__.py
│   ├── free_particle.py
│   ├── harmonic_oscillator.py        # HarmonicOscillator
│   └── potential_systems.py          # Puits, barrières
│
├── experiments/                      # Simulations complètes
│   ├── __init__.py
│   ├── base_experiment.py            # Classe Experiment abstraite
│   ├── wavepacket_evolution.py
│   └── measurement_statistics.py
│
├── validation/                       # Tests physiques
│   ├── __init__.py
│   ├── heisenberg_relations.py
│   ├── conservation_laws.py
│   └── ehrenfest_theorem.py
│
├── utils/                            # Outils transverses
│   ├── __init__.py
│   ├── numerical.py                  # FFT, différences finies, intégration
│   └── visualization.py              # Plots ψ, ρ, J
│
├── tests/                            # Tests unitaires (pytest)
│   ├── test_core/
│   ├── test_dynamics/
│   └── ...
│
└── examples/                         # Scripts exemples
    ├── example_gaussian_packet.py
    └── example_harmonic_oscillator.py
```


### 4.2 Flux de dépendances autorisées

```
experiments → systems → dynamics → core
    ↓           ↓          ↓         ↓
validation → peut accéder à toutes couches
    ↓
  utils (appelé par toutes couches)
```

**Règle stricte** : Pas de dépendance inverse (ex: `core` ne peut PAS importer `dynamics`)

### 4.3 Points d'entrée

1. **Mode script** : `examples/*.py` charge config, instancie Experiment, appelle `.run()`
2. **Mode interactif** : Import modules, construction manuelle objets
3. **Mode test** : `pytest tests/` vérifie invariants physiques

***

## 5. État actuel de l'implémentation

### 5.1 Implémenté (ou structure définie)

#### `core/constants.py`

- Classe `PhysicalConstants` pour charger h, ℏ, masses depuis yaml
- Méthode `validate_units()` pour vérifier ℏ = h/(2π)


#### `core/state.py`

- **Classe abstraite** `QuantumState` :
    - Méthodes abstraites : `norm()`, `normalize()`, `inner_product()`
    - Méthode concrète : `is_normalized()`
- **Classe concrète** `WaveFunctionState` :
    - Attributs : `spatial_grid`, `wavefunction` (np.ndarray complexe)
    - Méthodes : `probability_density()`, `probability_in_volume()`
    - **À implémenter** : calcul norme (intégrale discrète), produit scalaire
- **Classe concrète** `EigenStateBasis` :
    - Attributs : `eigenstates` (list), `coefficients`, `eigenvalues`
    - **À implémenter** : `validate_orthonormality()`, opérations algébriques


#### `core/operators.py`

- **Classe abstraite** `Observable` :
    - Méthodes abstraites : `apply()`, `expectation_value()`, `uncertainty()`, `eigensystem()`
    - Méthodes concrètes : `is_hermitian()`, `commutator()`
- **Classes concrètes partielles** :
    - `PositionOperator` : application = multiplication par r
    - `MomentumOperator` : application = -iℏ∇ (différences finies **à choisir**)
    - `Hamiltonian` : constructeur prend `mass`, `potential`, application = -ℏ²/2m Δ + V


#### `dynamics/evolution.py`

- **Classe** `TimeEvolution` :
    - Méthode `evolve_eigenstate()` : implémente Règle R3.3
    - Méthode `evolve_stationary_state()` : implémente Règle R3.4
    - Méthode `evolve_wavefunction()` : **STRUCTURE DÉFINIE**, schéma intégration **À CHOISIR**


#### `dynamics/measurement.py`

- **Classe** `QuantumMeasurement` :
    - Méthode `compute_probabilities()` : implémente Règle R2.2
    - Méthode `apply_reduction()` : implémente Règle R2.3
    - Méthode `measure_once()` : **À IMPLÉMENTER** (tirage aléatoire numpy)
    - Méthode `measure_ensemble()` : boucle sur `measure_once()`, **À IMPLÉMENTER**


#### `systems/harmonic_oscillator.py`

- **Classe** `HarmonicOscillator` :
    - Méthode `energy_eigenvalue(n)` : implémente Règle R6.1
    - Méthodes `creation_operator()`, `annihilation_operator()` : structure définie
    - **LIMITE** : Construction |0⟩ en représentation position nécessiterait fonction d'onde explicite (non dans extraits). **Solution adoptée** : travailler en base abstraite {|n⟩}


### 5.2 Prévu mais non codé / Récemment codé

#### Expériences

- `experiments/base_experiment.py` : classe `Experiment` abstraite avec méthodes `prepare_initial_state()`, `evolve_state()`, `perform_measurements()`, `run()`
- Implémentations concrètes : `WavePacketEvolution`, `MeasurementStatistics`


#### Validation

- `HeisenbergValidator.validate_position_momentum()` : calcule ΔX, ΔP, vérifie ≥ ℏ/2
- `ConservationValidator.validate_continuity_equation()` : calcule ∂ρ/∂t + ∇·J
- `EhrenfestValidator` : vérifie théorème R4.4


#### Utilitaires numériques

- `utils/numerical.py` : FFT, gradients, laplaciens, intégration
- `utils/visualization.py` : plots 1D/2D de ψ, ρ, J


### 5.3 Interfaces décidées (signatures clés)

```python
# core/state.py
class QuantumState:
    def norm(self) -> float: ...
    def inner_product(self, other: 'QuantumState') -> complex: ...
    def normalize(self) -> 'QuantumState': ...

# core/operators.py
class Observable:
    def apply(self, state: QuantumState) -> QuantumState: ...
    def expectation_value(self, state: QuantumState) -> float: ...
    def uncertainty(self, state: QuantumState) -> float: ...
    def eigensystem(self) -> tuple[np.ndarray, list[QuantumState]]: ...
    def commutator(self, other: 'Observable') -> 'Observable': ...

# dynamics/evolution.py
class TimeEvolution:
    def evolve_eigenstate(self, initial: EigenStateBasis, t0: float, t: float) -> EigenStateBasis: ...
    def evolve_wavefunction(self, initial: WaveFunctionState, t0: float, t: float, dt: float) -> WaveFunctionState: ...

# dynamics/measurement.py
class QuantumMeasurement:
    def compute_probabilities(self, state: QuantumState) -> dict[float, float]: ...
    def measure_once(self, state: QuantumState) -> tuple[float, QuantumState]: ...
    def apply_reduction(self, state: QuantumState, measured_value: float) -> QuantumState: ...

# experiments/base_experiment.py
class Experiment(ABC):
    def run(self) -> dict: ...
    def validate_physics(self) -> dict[str, bool]: ...
```


***

## 6. Gestion des expériences et simulations

### 6.1 Cycle type d'une expérience

```
┌────────────────────────────────────────┐
│ 1. PRÉPARATION (prepare_initial_state) │
│    - Créer |ψ(t₀)⟩                     │
│    - Vérifier normalisation            │
│    - Calculer propriétés initiales     │
└──────────────┬─────────────────────────┘
               ↓
┌────────────────────────────────────────┐
│ 2. DÉFINITION SYSTÈME                  │
│    (define_hamiltonian)                │
│    - Construire H                      │
│    - Vérifier hermiticité              │
│    - (Optionnel) Diagonaliser          │
└──────────────┬─────────────────────────┘
               ↓
┌────────────────────────────────────────┐
│ 3. ÉVOLUTION (evolve_state)            │
│    - Choisir méthode (R3.1/3.2/3.3)    │
│    - Intégrer de t₀ à t                │
│    - Vérifier conservation norme       │
└──────────────┬─────────────────────────┘
               ↓
┌────────────────────────────────────────┐
│ 4. MESURE (perform_measurements)       │
│    - Calculer ⟨A⟩, ΔA (R4.1, R4.2)     │
│    - Simuler mesures (R2.2, R2.3)      │
│    - Statistiques ensemble             │
└──────────────┬─────────────────────────┘
               ↓
┌────────────────────────────────────────┐
│ 5. ANALYSE (dans run ou post-process)  │
│    - Valider Heisenberg (R4.3)         │
│    - Valider Ehrenfest (R4.4)          │
│    - Valider conservation (R5.1, R5.2) │
└────────────────────────────────────────┘
```


### 6.2 Contraintes par étape

**Étape 1 (Préparation)** :

- Autorisé : Construire `WaveFunctionState` ou `EigenStateBasis`
- Obligatoire : Vérifier `⟨ψ|ψ⟩ = 1`
- Source yaml : Paramètres états (largeur gaussienne, position initiale, etc.)

**Étape 2 (Hamiltonien)** :

- Autorisé : Instancier `Hamiltonian` avec masse, potentiel
- Obligatoire : Vérifier `H† = H`
- Source yaml : Paramètres système (m, ω, forme potentiel)

**Étape 3 (Évolution)** :

- Autorisé : Appeler `TimeEvolution.evolve_*()` selon type état
- Interdit : Modifier H en cours d'évolution (sauf H(t) explicite)
- Vérification continue : `⟨ψ(t)|ψ(t)⟩ = 1` à chaque pas

**Étape 4 (Mesure)** :

- Autorisé : `QuantumMeasurement` sur état évolué
- Pour statistiques : Répéter sur N copies état initial (pas sur état post-mesure)
- Respecter réduction paquet d'ondes (R2.3)

**Étape 5 (Analyse)** :

- Calculer indicateurs physiques
- Comparer à prédictions théoriques (Heisenberg, Ehrenfest)
- Générer rapport validation


### 6.3 Exemple concret : `WavePacketEvolution`

**Objectif** : Observer étalement paquet gaussien libre, vérifier Heisenberg

**Étapes** :

1. Préparation : Gaussienne centrée x₀, largeur σₓ, impulsion k₀
2. Hamiltonien : H = P²/2m (V=0, particule libre)
3. Évolution : De t=0 à t=T par pas dt
4. Mesures : À chaque temps, calculer ⟨X⟩, ⟨P⟩, ΔX, ΔP
5. Analyse :
    - Vérifier ΔX(t)·ΔP(t) ≥ ℏ/2 pour tout t
    - Vérifier d⟨X⟩/dt = ⟨P⟩/m (Ehrenfest)
    - Observer ΔX(t) croît avec t (étalement)

**Paramètres yaml** :

```yaml
experiments:
  wavepacket_evolution:
    initial_state:
      type: "gaussian"
      x0: 0.0
      sigma_x: 1.0e-9  # m
      k0: 1.0e10       # m^-1
    system:
      mass: 9.11e-31   # kg (électron)
      potential: "free"
    evolution:
      t_initial: 0.0
      t_final: 1.0e-15  # s
      dt: 1.0e-17
      times_sample: [0, 5e-16, 1e-15]  # temps pour mesures
```


***

## 7. Configuration et paramètres (parameters.yaml)

### 7.1 Rôle du fichier

**Objectif** : Centraliser TOUS paramètres physiques et numériques ajustables
**Principe** : Code = logique, YAML = valeurs
**Chargement** : Au démarrage expérience, une seule fois
**Validation** : Vérifier cohérence (ex: h/(2π) = ℏ) à chargement

### 7.2 Structure décidée

```yaml
# ============================================================
# CONSTANTES PHYSIQUES FONDAMENTALES
# ============================================================
physical_constants:
  h: 6.62607015e-34        # Constante Planck (J·s)
  hbar: 1.054571817e-34    # ℏ = h/(2π) (J·s)
  c: 2.99792458e8          # Vitesse lumière (m/s)
  m_electron: 9.1093837015e-31   # Masse électron (kg)
  m_proton: 1.67262192369e-27    # Masse proton (kg)
  e_charge: 1.602176634e-19      # Charge élémentaire (C)
  epsilon_0: 8.8541878128e-12    # Permittivité vide (F/m)

# ============================================================
# PARAMÈTRES NUMÉRIQUES GÉNÉRAUX
# ============================================================
numerical_parameters:
  spatial_discretization:
    dimension: 1             # 1, 2 ou 3
    grid_type: "uniform"     # Pour l'instant, uniquement uniforme
    nx: 1024                 # Nombre points direction x
    x_min: -1.0e-8           # Borne inférieure (m)
    x_max: 1.0e-8            # Borne supérieure (m)
    # ny, nz, y_min, ... si dimension > 1
    
  temporal_discretization:
    dt: 1.0e-17              # Pas temps (s)
    t_initial: 0.0
    t_final: 1.0e-15
    
  tolerances:
    normalization_check: 1.0e-10
    hermiticity_check: 1.0e-10
    orthonormality_check: 1.0e-8
    heisenberg_inequality: 1.0e-10
    conservation_probability: 1.0e-9

# ============================================================
# PARAMÈTRES SYSTÈMES PHYSIQUES
# ============================================================
systems:
  free_particle:
    mass: 9.1093837015e-31   # m électron (kg)
    
  harmonic_oscillator:
    mass: 9.1093837015e-31
    omega: 1.0e15            # Pulsation (rad/s)
    n_max_fock: 50           # Troncature base {|n⟩}
    
  potential_systems:
    # Configuration spécifique selon type
    infinite_well:
      width: 1.0e-9          # Largeur puits (m)
    # Autres potentiels...

# ============================================================
# CONFIGURATION EXPÉRIENCES
# ============================================================
experiments:
  wavepacket_evolution:
    initial_state:
      type: "gaussian"
      x0: 0.0
      sigma_x: 1.0e-9
      k0: 1.0e10
    evolution_times: [0, 5e-16, 1e-15]  # Temps échantillonnage
    
  measurement_statistics:
    observable: "position"    # "position", "momentum", "energy"
    n_measurements: 10000
    initial_state:
      type: "custom"
      # Définition spécifique...

# ============================================================
# VISUALISATION (optionnel)
# ============================================================
visualization:
  plot_wavefunction: true
  plot_probability_density: true
  plot_current_density: false
  save_figures: true
  output_directory: "./results/"
  dpi: 150
```


### 7.3 Catégories obligatoires vs optionnelles

**OBLIGATOIRE** (expérience ne peut se lancer sans) :

- `physical_constants.h`, `hbar` (cohérence validée)
- `numerical_parameters.tolerances.*`
- `systems.<type_système>` pour système utilisé
- `experiments.<nom_expérience>` pour expérience lancée

**OPTIONNEL** :

- Constantes non utilisées (ex: `c` si pas relativiste)
- Paramètres visualisation
- Paramètres systèmes non utilisés

**INTERDIT dans YAML** (doit rester dans code) :

- Choix algorithmes (schéma intégration, méthode FFT)
- Structure classes
- Logique métier

***

## 8. Limites actuelles et points ouverts

### 8.1 Limites théoriques (imposées par le cours)

#### L1 : Pas de méthode numérique d'intégration fournie

**Problème** : Règle R3.2 donne `iℏ∂ψ/∂t = Hψ` mais aucun schéma résolution
**Impact** : `TimeEvolution.evolve_wavefunction()` nécessite choix algorithme externe au cours
**Options** : Euler explicite (instable), Crank-Nicolson (stable), split-operator (efficace FFT)
**Décision requise** : Implémenter un schéma, justifier stabilité/précision

#### L2 : Fonctions d'onde explicites partielles

**Problème** : États fondamentaux non donnés en représentation position
**Exemples** :

- Oscillateur harmonique |0⟩ : Complément BV mentionné mais extraits absents
- Atome H : Seulement E₀, a₀ fournis, pas fonctions radiales complètes
**Impact** :
- HO : **Solution adoptée** = travailler en base abstraite {|n⟩}
- Atome H : **Non implémentable** sans compléments théoriques


#### L3 : Traitement spectre continu incomplet

**Problème** : Règle R2.2 pour spectre continu donnée formellement, mais pas discrétisation
**Impact** : Mesure position/impulsion (spectre continu) nécessite approximation par grille
**Solution partielle** : Grille = spectre discret effectif, justifier convergence

#### L4 : Pas de théorie spin dans extraits

**Problème** : Chapitre IV (spin 1/2) mentionné mais non fourni
**Impact** : **Impossible** d'implémenter systèmes avec spin (électron réel, Stern-Gerlach)

#### L5 : Pas de particules identiques

**Problème** : Tome I fait "première incursion" mais détails dans Tome II
**Impact** : Pas de systèmes multi-particules avec statistiques quantiques

### 8.2 Limites numériques actuelles

#### N1 : Discrétisation spatiale

**Choix actuel** : Grille uniforme, différences finies pour dérivées
**Alternative non implémentée** : FFT (efficace mais conditions périodiques), éléments finis
**Paramètres** : nx, xmin, xmax (yaml)
**Test requis** : Convergence quand nx→∞

#### N2 : Représentation opérateurs matriciels

**Problème** : Impulsion P en représentation position = dérivée, pas matrice directe
**Solution actuelle** : Application = différences finies sur tableau `wavefunction`
**Alternative** : FFT (P diagonal en représentation impulsion)

#### N3 : Calcul valeurs/vecteurs propres

**Méthode** : `np.linalg.eigh()` pour matrices hermitiennes
**Limite** : Dimension max ~10⁴×10⁴ (RAM)
**Oscillateur** : Matrice creuse (tri-diagonale), méthodes spécialisées possibles

#### N4 : Intégration numérique (probabilités)

**Problème** : `∫|ψ|² d³r` approximé par somme discrète
**Formule** : `Σᵢ |ψᵢ|² · dV` où dV = dx·dy·dz
**Erreur** : O(dx²) si différences finies ordre 2

### 8.3 Points ouverts nécessitant décision

#### D1 : Schéma intégration temporelle

**Question** : Quel algorithme pour `evolve_wavefunction()` ?
**Options** :

1. **Euler explicite** : `ψ(t+dt) = ψ(t) - i(dt/ℏ)Hψ(t)`
    - Avantage : Simple
    - Inconvénient : Instable, ne conserve pas norme
2. **Crank-Nicolson** : `(1+iHdt/2ℏ)ψ(t+dt) = (1-iHdt/2ℏ)ψ(t)`
    - Avantage : Stable, conserve norme
    - Inconvénient : Résolution système linéaire chaque pas
3. **Split-operator** : `ψ(t+dt) = exp(-iV·dt/ℏ) · FFT⁻¹[exp(-iP²dt/2mℏ) · FFT[ψ(t)]]`
    - Avantage : Rapide (FFT), conserve norme
    - Inconvénient : Conditions périodiques implicites

**Recommandation initiale** : Crank-Nicolson (stabilité prioritaire), optimiser plus tard

#### D2 : Calcul gradient/laplacien

**Question** : Ordre différences finies ?
**Options** :

- Ordre 2 : `∂ψ/∂x ≈ (ψᵢ₊₁ - ψᵢ₋₁)/(2dx)`
- Ordre 4 : Plus précis, plus coûteux
**Décision provisoire** : Ordre 2 par défaut, configurable si besoin


#### D3 : Gestion bords grille spatiale

**Problème** : ψ défini en x_min, x_max, que faire aux bords ?
**Options** :

1. Conditions Dirichlet : ψ(x_min) = ψ(x_max) = 0 (puits infini implicite)
2. Conditions périodiques : ψ(x_min) = ψ(x_max) (pour FFT)
3. Conditions absorbantes : éviter réflexions
**Décision provisoire** : Dirichlet (plus simple), documenter dans docstring

#### D4 : Construction état fondamental oscillateur

**Problème** : a|0⟩ = 0 définit |0⟩ abstraitement, comment obtenir ψ₀(x) ?
**Option 1** : Utiliser formule analytique (polynômes Hermite) - **hors extraits cours**
**Option 2** : Résoudre numériquement équation aux valeurs propres pour n=0
**Option 3** : Travailler uniquement en base {|n⟩} abstraite (matrices)
**Décision adoptée** : **Option 3** (cohérent avec cours fourni)

#### D5 : Tirage aléatoire mesures

**Implémentation** : `np.random.choice(eigenvalues, p=probabilities)`
**Graine** : Paramètre `random_seed` optionnel pour reproductibilité
**Question ouverte** : Faut-il logger séquence mesures ou seulement statistiques finales ?

### 8.4 Extensions futures (hors périmètre actuel)

#### E1 : Systèmes 2D/3D

**Requis** :

- Grilles 2D/3D (déjà prévu dans yaml)
- Laplacien 2D/3D
- Visualisations adaptées (contours, isosurfaces)
**Effort** : Modéré, généralisation code 1D


#### E2 : Potentiels dépendant du temps V(r,t)

**Requis** :

- Modifier `Hamiltonian.__init__()` pour accepter V(r,t)
- Adapter `TimeEvolution` (pas de décomposition spectrale simple)
**Effort** : Faible si intégration numérique déjà implémentée


#### E3 : Systèmes multi-particules (sans spin/identité)

**Exemple** : 2 particules sans interaction → produit tensoriel états
**Requis** : Module `core/hilbert_space.py` pour produits tensoriels
**Effort** : Moyen, complexité augmente vite

#### E4 : Atome hydrogène complet

**Bloqueurs actuels** :

- Résolution équation radiale non fournie
- Fonctions Laguerre, harmoniques sphériques (Complément AVI mentionné)
**Requis** : Accès Compléments ou implémentation externe fonctions spéciales
**Effort** : Important (nécessite théorie additionnelle)


#### E5 : Spin et systèmes 2 niveaux

**Bloqueurs** : Chapitre IV non fourni dans extraits
**Requis** : Contenu Chapitre IV complet
**Effort** : Moyen une fois théorie disponible

***

## 9. Références traçabilité

### 9.1 Tableau correspondance Règles ↔ Sources cours

| Règle | Description courte | Source cours |
| :-- | :-- | :-- |
| R1.1 | Planck-Einstein | [file:1, Chap I, §A-1, équations A-1,A-2] |
| R1.2 | De Broglie | [file:1, Chap I, §B-1] |
| R1.3 | Commutateurs canoniques | [file:1, Chap III, §B-5-a, équations B-33] |
| R2.1 | Densité probabilité | [file:1, Chap I, §B-2] |
| R2.2 | Probabilités mesure | [file:1, Chap III, §B-3-b, équations B-4,B-7] |
| R2.3 | Réduction paquet | [file:1, Chap III, §B-3-c, équations B-30,B-31] |
| R3.1 | Schrödinger abstrait | [file:1, Chap III, §B-4] |
| R3.2 | Schrödinger position | [file:1, Chap I, §B-2, équation D-1] |
| R3.3 | Décomposition spectrale | [file:1, Chap III, §D-2-a, équations D-54] |
| R3.4 | États stationnaires | [file:1, Chap III, §D-2-b, équation D-57] |
| R4.1 | Valeur moyenne | [file:1, Chap III, §C-4] |
| R4.2 | Écart quadratique | [file:1, Chap III, §C-5] |
| R4.3 | Heisenberg | [file:1, Chap III, §C-5, inégalités] |
| R4.4 | Ehrenfest | [file:1, Chap III, §D-1-d] |
| R4.5 | Hermiticité | [file:1, Chap II, §D-1] |
| R5.1 | Conservation norme | [file:1, Chap III, §D-1-c] |
| R5.2 | Équation continuité | [file:1, Chap III, §D-1-c] |
| R6.1 | Spectre HO | [file:1, Chap V, §B] |
| R6.2 | Algèbre a,a† | [file:1, Chap V, §B] |
| R6.3 | Action échelle | [file:1, Chap V, §C] |

### 9.2 Tableau correspondance Classes ↔ Règles

| Classe | Règles implémentées | Tests requis |
| :-- | :-- | :-- |
| `QuantumState` | - | Normalisation |
| `WaveFunctionState` | R2.1 | ∫\|ψ\|² = 1 |
| `EigenStateBasis` | R2.2 | Orthonormalité base |
| `Observable` | R4.1, R4.2, R4.5 | Hermiticité |
| `PositionOperator` | R1.3 (avec P) | [X,P] = iℏ |
| `MomentumOperator` | R1.3 (avec X) | P = -iℏ∇ |
| `Hamiltonian` | R3.1, R3.2 | H† = H |
| `TimeEvolution` | R3.1, R3.3, R3.4, R5.1 | ⟨ψ\|ψ⟩ constant |
| `QuantumMeasurement` | R2.2, R2.3 | Σ P(aₙ) = 1 |
| `HarmonicOscillator` | R6.1, R6.2, R6.3 | [a,a†]=1, Eₙ correct |
| `HeisenbergValidator` | R4.3 | ΔX·ΔP ≥ ℏ/2 |
| `ConservationValidator` | R5.2 | ∂ρ/∂t + ∇·J = 0 |
| `EhrenfestValidator` | R4.4 | d⟨R⟩/dt = ⟨P⟩/m |


***

## 10. Checklist implémentation future

### Pour chaque nouvelle classe/méthode

- [ ] **Traçabilité** : Docstring mentionne règle(s) implémentée(s) et source cours
- [ ] **Invariants** : Tests unitaires vérifient propriétés physiques (hermiticité, normalisation, etc.)
- [ ] **Paramètres** : Valeurs numériques viennent de yaml, pas en dur
- [ ] **Exceptions** : Lever erreur explicite si préconditions non respectées (ex: état non normé)
- [ ] **Documentation** : Hypothèses numériques explicitées (ordre différences finies, etc.)


### Pour chaque expérience

- [ ] **Config yaml** : Section dédiée avec tous paramètres
- [ ] **Cycle complet** : Préparation → Évolution → Mesure → Analyse
- [ ] **Validation** : Appel méthodes validation (Heisenberg, conservation, Ehrenfest)
- [ ] **Résultats** : Export structuré (dict ou fichier) avec métadonnées (params, temps calcul)
- [ ] **Visualisation** : Au moins un plot résumant résultats physiques


### Tests physiques obligatoires

1. **Test normalisation** : Tous états manipulés ont ⟨ψ|ψ⟩ ≈ 1
2. **Test hermiticité** : Toutes observables vérifient A† = A
3. **Test Heisenberg** : ΔX·ΔP ≥ ℏ/2 pour états testés
4. **Test conservation** : Norme constante durant évolution
5. **Test commutateurs** : [X,P] = iℏ (sur états tests)

***

## 11. Glossaire des symboles

| Symbole | Signification | Unité SI |
| :-- | :-- | :-- |
| h | Constante Planck | J·s |
| ℏ | Constante réduite (h/2π) | J·s |
| m | Masse | kg |
| E | Énergie | J |
| p | Impulsion | kg·m/s |
| k | Vecteur d'onde | m⁻¹ |
| ω | Pulsation | rad/s |
| λ | Longueur d'onde | m |
| ν | Fréquence | Hz |
| ψ(r,t) | Fonction d'onde | m⁻³/² (3D) |
| ρ(r,t) | Densité probabilité | m⁻³ (3D) |
| J(r,t) | Courant probabilité | m⁻²·s⁻¹ (3D) |
| \|ψ⟩ | Vecteur d'état (ket) | Sans dimension (normé) |
| ⟨ψ\| | Vecteur dual (bra) | Sans dimension |
| A | Observable générique | Variable (dépend grandeur) |
| H | Hamiltonien | J |
| R, X, Y, Z | Opérateurs position | m |
| P, Pₓ, Pᵧ, Pᵤ | Opérateurs impulsion | kg·m/s |
| a, a† | Annihilation, création (HO) | Sans dimension |
| n | Nombre quantique | Entier ≥ 0 |
| ΔA | Écart quadratique moyen A | Unité de A |
| [A,B] | Commutateur | Unité de AB |


***

## 12. Synthèse finale pour l'agent de code

### Ce que l'agent DOIT respecter absolument

1. **Toute équation implémentée** provient d'une règle R*.* traçable au cours
2. **Aucune extrapolation** : si formule manquante, demander clarification
3. **Structure modulaire** : respecter architecture dossiers et flux dépendances
4. **Paramètres externes** : tous ajustables via `parameters.yaml`
5. **Tests physiques** : chaque classe/méthode doit avoir tests unitaires validant invariants physiques
6. **Documentation complète** : docstrings avec références règles et sources
