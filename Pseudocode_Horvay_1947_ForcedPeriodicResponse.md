# Pseudocode — Horvay (1947) "Rotor Blade Flapping Motion"
## Forced Periodic Response Only

**Source:** Horvay, G. (1947). *Rotor Blade Flapping Motion.*
Quarterly of Applied Mathematics, Vol. V, No. 2, pp. 149–167.

**Scope of this pseudocode:**
This pseudocode covers the forced periodic (steady-state) flapping
response only, following Horvay's original derivation exactly.
The governing equation is written in the **time domain** as in the
original paper. All equation numbers refer to Horvay (1947).

**Fidelity level:** Rigid blade, 1-DOF flap, linear quasi-steady
strip-theory aerodynamics, uniform inflow, centrally hinged
(hinge at rotor axis), small-angle assumption throughout.

---

## 1. Symbol Legend

All symbols follow Horvay (1947) exactly.

### 1.1 Independent Variable

| Symbol | Meaning | Units |
|---|---|---|
| t | Time | s |
| ψ | Blade azimuth angle, ψ = ωt (measured from rearmost position) | rad |

### 1.2 Dependent Variable

| Symbol | Meaning | Units |
|---|---|---|
| β(t) | Flapping angle — angle between blade and rotor plane, positive upward | rad |
| β̇ | dβ/dt — flapping angular rate | rad/s |
| β̈ | d²β/dt² — flapping angular acceleration | rad/s² |

### 1.3 Primary Physical Parameters (Direct Inputs)

| Symbol | Meaning | Units |
|---|---|---|
| ρ | Air density | kg/m³ |
| c | Blade chord (constant along span) | m |
| a | Lift-curve slope of blade section | 1/rad |
| R | Blade length (radius) | m |
| I | Flapping moment of inertia of blade about hinge, I = ∫₀ᴿ r² dm | kg·m² |
| ω | Rotor angular velocity | rad/s |
| m | Blade mass | kg |
| g | Gravitational acceleration | m/s² |
| r_cg | Distance of blade centre of gravity from flapping hinge | m |

### 1.4 Operating Condition Parameters (Direct Inputs)

| Symbol | Meaning | Units |
|---|---|---|
| μ | Advance ratio — ratio of forward flight speed to blade tip speed | — |
| λ | Inflow ratio — (sinking speed − induced velocity) / (ωR) | — |

### 1.5 Control Inputs (Constant at a Trim Condition)

| Symbol | Meaning | Units |
|---|---|---|
| ϑ₀ | Collective pitch angle | rad |
| ϑ_c | Lateral cyclic pitch amplitude | rad |
| ϑ_s | Fore-and-aft cyclic pitch amplitude (−ϑ_s is fore-and-aft setting) | rad |
| ϑ(t) | Total blade pitch: ϑ = ϑ₀ + ϑ_c cos ωt + ϑ_s sin ωt | rad |

### 1.6 Derived Non-Dimensional Groups

| Symbol | Meaning | Definition | Units |
|---|---|---|---|
| n | Aerodynamic damping coefficient | n = ρcaR⁴ / (8I) | — |

### 1.7 Velocity Components at Blade Element (Eq. 3a, 3b)

| Symbol | Meaning | Units |
|---|---|---|
| r | Spanwise position of blade element from hinge | m |
| U_T | Tangential component of air velocity at element | m/s |
| U_P | Perpendicular component of air velocity at element | m/s |
| φ | Inflow angle at blade element, φ = U_P / U_T (small angle) | rad |

### 1.8 Governing Equation Coefficients (Eq. 5a, 5b, 5c)

| Symbol | Meaning | Units |
|---|---|---|
| p(t) | Time-varying aerodynamic damping coefficient | 1/s |
| s(t) | Time-varying aerodynamic spring (stiffness) coefficient | 1/s² |
| E(ωt) | Aerodynamic and gravitational forcing function | rad/s² |

### 1.9 Forced Response Fourier Coefficients (Eq. 10)

| Symbol | Meaning | Units |
|---|---|---|
| a₀ | Coning angle — constant (mean) flapping angle | rad |
| a₁ | Amplitude of 1/rev cosine flapping term | rad |
| b₁ | Amplitude of 1/rev sine flapping term | rad |
| a₂ | Amplitude of 2/rev cosine flapping term | rad |
| b₂ | Amplitude of 2/rev sine flapping term | rad |
| aₙ, bₙ | Amplitudes of n/rev cosine and sine flapping terms | rad |

### 1.10 Numerical Control

| Symbol | Meaning |
|---|---|
| N | Truncation order — number of harmonics retained in the Fourier expansion |
| tol | Convergence tolerance for iterative harmonic refinement |

---

## 2. Parameter Extraction Block

```
PROCEDURE setup_horvay_problem(inputs):

    # ----------------------------------------------------------------
    # 2.1  Read all direct physical inputs
    # ----------------------------------------------------------------
    READ ρ              # air density [kg/m³]
    READ c              # blade chord [m]
    READ a              # lift-curve slope [1/rad]
    READ R              # blade radius [m]
    READ I              # flapping moment of inertia [kg·m²]
    READ ω              # rotor angular velocity [rad/s]
    READ m              # blade mass [kg]
    READ g              # gravitational acceleration [m/s²]
    READ r_cg           # c.g. distance from hinge [m]

    # ----------------------------------------------------------------
    # 2.2  Read operating condition (direct inputs, not derived)
    # ----------------------------------------------------------------
    READ μ              # advance ratio [-]
    READ λ              # inflow ratio [-]

    # ----------------------------------------------------------------
    # 2.3  Read control settings (constant at trim point)
    # ----------------------------------------------------------------
    READ ϑ₀             # collective pitch [rad]
    READ ϑ_c            # lateral cyclic pitch [rad]
    READ ϑ_s            # fore-and-aft cyclic pitch [rad]

    # ----------------------------------------------------------------
    # 2.4  Compute derived non-dimensional aerodynamic parameter
    # ----------------------------------------------------------------
    n ← ρ * c * a * R⁴ / (8 * I)          # aerodynamic damping coeff. (Eq. 7)

    # ----------------------------------------------------------------
    # 2.5  Compute blade weight moment term for forcing function
    # ----------------------------------------------------------------
    weight_moment ← m * g * r_cg / I       # appears in E(ωt), Eq. 5c

    # ----------------------------------------------------------------
    # 2.6  Read numerical control parameters
    # ----------------------------------------------------------------
    READ N              # harmonic truncation order (typical: 3–5 for μ ≤ 0.35)
    READ tol            # convergence tolerance (e.g. 1e-6)

    # ----------------------------------------------------------------
    # 2.7  Assemble and return parameter struct
    # ----------------------------------------------------------------
    RETURN P = {
        n,              # aerodynamic damping coefficient
        ω,              # rotor angular velocity
        μ,              # advance ratio
        λ,              # inflow ratio
        ϑ₀,             # collective pitch
        ϑ_c,            # lateral cyclic
        ϑ_s,            # fore-and-aft cyclic
        weight_moment,  # mgr_cg/I
        N,              # harmonic truncation order
        tol             # convergence tolerance
    }
END
```

---

## 3. Governing Equation — Hill's Differential Equation in Time Domain

### 3.1 Full Forced Equation (Eq. 5, 5a)

The flapping motion of the blade is governed by (Eq. 5):

```
J(β) = E(ωt)
```

where the operator J is defined as (Eq. 5a):

```
J(β) ≡ β̈  + p(t)·β̇ + s(t)·β
```

### 3.2 Time-Varying Coefficients (Eq. 5b)

The aerodynamic damping coefficient p(t) and spring coefficient s(t) are:

```
p(t) = n·ω · [1 + (4/3)·μ·sin(ωt)]

s(t) = ω² · [1 + (4/3)·n·μ·cos(ωt) + n·μ²·sin(2ωt)]
```

### 3.3 Forcing Function (Eq. 5c)

The aerodynamic and gravitational forcing is:

```
E(ωt) = − weight_moment
         + n·ω² · {
               (4/3)·λ  +  (1 + μ²)·ϑ
             + [2·λ·μ  +  (8/3)·μ·ϑ] · sin(ωt)
             −  μ²·ϑ · cos(2ωt)
           }
```

where `weight_moment = mgr_cg/I`.

**Structure of E(ωt):**
- The first term is the constant moment of the blade weight about the
  hinge (Horvay notes this is usually negligibly small)
- The remaining terms are the aerodynamic excitation, comprising:
  - A constant (mean) part driven by λ and ϑ₀
  - A 1/rev sine term driven by forward speed μ, inflow λ, and
    total cyclic pitch ϑ
  - A 2/rev cosine term driven by μ² and total cyclic pitch ϑ

### 3.4 Blade Pitch Variation (Eq. 5d)

The blade pitch angle varies with azimuth as:

```
ϑ = ϑ(t) = ϑ₀ + ϑ_c·cos(ωt) + ϑ_s·sin(ωt)
```

where ϑ₀ is the collective pitch, ϑ_c the lateral, and −ϑ_s the
fore-and-aft control setting.

### 3.5 Velocity Components (Eq. 3a, 3b)

For reference, the air velocity components at blade element dm at
spanwise position r are:

```
U_T = r·ω + μ·ω·R·sin(ωt)          # tangential component

U_P = λ·ω·R − r·β̇ − μ·ω·R·β·cos(ωt)   # perpendicular component
```

### 3.6 Homogeneous Equation

Setting E(ωt) = 0 gives the homogeneous (unforced) equation:

```
J(β) = 0
```

i.e.:

```
β̈ + p(t)·β̇ + s(t)·β = 0
```

This governs the free (transient) flapping motion.

### 3.7 General Solution Structure (Eq. 8)

The general solution of the full equation is:

```
β(t) = β₀(t) + h₁·β₁(t) + h₂·β₂(t)
```

where:
- β₀(t) is the **particular integral** — the forced periodic
  (steady-state) response
- β₁(t), β₂(t) are the two independent solutions of the homogeneous
  equation — the transient modes
- h₁, h₂ are scalar coefficients fixed by initial conditions
  β(0) and β_dot(0)

---

## 4. Forced Periodic Response β₀(t)

### 4.1 Fourier Ansatz (Eq. 10)

Because E(ωt) is periodic with period 2π/ω, the steady-state
response β₀(t) is also periodic with the same period. It is assumed
in the Fourier series form:

```
β₀(t) = a₀
       + a₁·cos(ωt)  + b₁·sin(ωt)
       + a₂·cos(2ωt) + b₂·sin(2ωt)
       + a₃·cos(3ωt) + b₃·sin(3ωt)
       + ...
```

**Physical interpretation of the leading coefficients:**
- a₀ : coning angle — mean out-of-plane blade position
- a₁ : longitudinal tilt of the rotor cone axis (forward when ω
       is counterclockwise)
- b₁ : lateral tilt of the rotor cone axis (to port side when ω
       is counterclockwise)
- a₂, b₂, a₃, b₃, ... : higher harmonics — describe blade motion
       in and out of the cone surface

### 4.2 Linear System of Equations for the Fourier Coefficients (Eq. 11)

Substituting the ansatz (Eq. 10) into the full equation J(β) = E(ωt)
and equating Fourier coefficients of like harmonics (1, cos ψ, sin ψ,
cos 2ψ, sin 2ψ, ...) yields the infinite linear system (Eq. 11):

```
Constant (1):
    a₀  +  (1/2)·μ²·n·b₂
        = − weight_moment/ω²
          + n·{ ½·λ  +  (1 + μ²)·ϑ₀ }

cos(ωt):
    (4/3)·μ·n·a₀  +  (1 + ½·μ²)·b₁  −  (2/3)·μ·a₂  +  (1/2)·n·μ²·b₃
        = (1 + (1/2)·μ²)·ϑ_c    

sin(ωt):
    −(1 − (1/2)μ²)·a₁  −  (2/3)·μ·b₂  −  (1/2)·μ²·a₃
        = 2·λ·μ  +  (8/3)·μ·ϑ₀  +  (1 + 3/2·μ²)·ϑ_s

cos(2ωt):
    (4/3)·μ·n·a₁  −  3·a₂  +  2·n·b₂  −  (4/3)·n·μ·a₃  + ...
        = − μ²·n·ϑ₀  −  (4/3)·μ·n·ϑ_s

sin(2ωt):
    μ²·n·a₀  +  (4/3)·μ·n·b₁  −  2·n·a₂  −  3·b₂  −  (4/3)·μ·n·b₃  + ...
        = (4/3)·μ·n·ϑ_c

...and so on for cos(3ωt), sin(3ωt), etc.
```

**Note:** The full infinite system is truncated by setting
aₙ = bₙ = 0 for n > N. The first (2N + 1) equations are then
solved for a₀, a₁, b₁, ..., a_N, b_N.

### 4.3 Solution Procedure

```
FUNCTION solve_forced_response(P):

    # 1. Set truncation order
    N ← P.N (P is the parameter struct from section ##2)

    # 2. Build unknown vector
    x = [a₀, a₁, b₁, a₂, b₂, ..., a_N, b_N]    # length = 2N + 1

    # 3. Assemble the linear system A·x = b from Eq. 11
    #    by substituting the Fourier ansatz into J(β) = E(ωt)
    #    and collecting coefficients of each harmonic

    FOR each harmonic k = 0, 1, ..., N:
        FOR each unknown harmonic j = 0, 1, ..., N:
            A[k, j] ← coefficient of (aⱼ or bⱼ) in the
                       k-th harmonic equation of Eq. 11
        b[k] ← right-hand side of the kth harmonic equation
                from E(ωt), using P.n, P.μ, P.λ,
                P.ϑ₀, P.ϑ_c, P.ϑ_s, P.weight_moment, P.ω
    END FOR

    # 4. Solve the linear system
    x ← LINEAR_SOLVE(A, b)          # standard Gaussian elimination

    # 5. Extract Fourier coefficients
    a₀ ← x[0]
    a₁ ← x[1];    b₁ ← x[2]
    a₂ ← x[3];    b₂ ← x[4]
    ...
    a_N ← x[2N-1]; b_N ← x[2N]

    # 6. Optional: increase N until convergence
    REPEAT
        N ← N + 1
        recompute A, b, x
    UNTIL max(|Δa₀|, |Δa₁|, |Δb₁|) < P.tol

    # 7. Reconstruct β₀(t) on a fine time grid for output and plotting
    FOR each time sample tₖ in [0, 2π/ω):
        β₀(tₖ) ← a₀
                  + Σₙ [ aₙ·cos(n·ω·tₖ) + bₙ·sin(n·ω·tₖ) ]
    END FOR

    RETURN {
        a₀,                  # coning angle
        a₁, b₁,              # 1/rev flapping coefficients
        a₂, b₂,              # 2/rev flapping coefficients
        higher harmonics,    # aₙ, bₙ for n > 2
        β₀(t)                # reconstructed time history
    }
END
```

### 4.4 Verification Example from Horvay (Eq. 12a, 12b)

Horvay provides a worked numerical example with the following
input parameters (Eq. 12a):

```
n              = 1.7
mgr_cg / (I·ω²) = 0.03          # normalized weight moment
λ              = −0.10
ϑ₀             = 0.2  rad
ϑ_c            = 0.0  rad        # no lateral cyclic
ϑ_s            = 0.0  rad        # no fore-and-aft cyclic
μ              = 0.34738
```

The resulting forced (steady-state) flapping response is (Eq. 12b):

```
β₀(ψ) = + 0.124
         − 0.125·cos(ψ)
         − 0.057·sin(ψ)
         − 0.012·cos(2ψ)
         + 0.007·sin(2ψ)
         − 0.001·cos(3ψ)
         + ...
```

**Interpretation of this result:**
- a₀ = +0.124 rad : the blade cones upward at ~7.1°
- a₁ = −0.125 rad : longitudinal tilt of the disc (~7.2° amplitude)
- b₁ = −0.057 rad : lateral tilt of the disc (~3.3° amplitude)
- a₂ = −0.012 rad : 2/rev cosine harmonic (small)
- b₂ = +0.007 rad : 2/rev sine harmonic (small)
- a₃ = −0.001 rad : 3/rev harmonic (negligible)

**Convergence note (from Horvay Section 2):**
For μ ≤ 0.35, termination of the series at the second harmonic
(N = 2) gives sufficiently accurate results. The convergence of
the series was confirmed by Lock, Wheatley and others from
numerous numerical examples.

**Implementation check:**
When implementing the solution procedure in Section 4.3, the
output Fourier coefficients must match the values in Eq. 12b
to within the convergence tolerance tol. If they do not match,
check:
1. The sign convention of ϑ in Eq. 5d (positive signs on both
   cyclic terms)
2. The coefficient of sin(ωt) in p(t) — should be 4/3·μ
3. The inclusion of the weight_moment term in the constant equation
4. The 2/rev cosine term in E(ωt) — driven by μ²·ϑ
```

---

*End of pseudocode — Horvay (1947) forced periodic response.*

*Stability analysis (transient solutions β₁, β₂ and apparent*
*damping n_app) is covered separately following Horvay Sections 3–5.*
