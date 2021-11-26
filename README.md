# tumour-automata
**A cellular automata model of tumour cell proliferation**

Modelling cell population dynamics can give useful insight into behaviour in a physiological context
Cellular automata models analyse the evolution of cells on a matrix given a set of starting conditions and rules for cell behaviour over time
In this project, a simple 2D cellular automaton model was developed for cancer cell proliferation that evolves based on known starting conditions, coefficients of motility, reproduction and cell death
The baseline forward model was expanded to incorporate:
Immune cell agents; and
Nutrient diffusion 
Sequential Monte Carlo Approximate Bayesian Computation (SMC-ABC) was subsequently applied to infer parameter values from blinded traces

**Pseudocode:**

Seed a 50x50 grid with 10 cancer cells
Optionally: add immune cells at n=agent_ratio*n_seed (default=2)
For n_iterations (default 3), while t<40 timesteps:
  Find non-zero cell coordinates
    For each cancer cell:
    Find nonzero neighbours and gaps
    If n_gaps<=4, cell dies with probability ∝p(mortality)
    If n_gaps>0, cell moves with probability∝p(multiply) and p(motility)
  (For immune cells):
    Find nonzero neighbours and gaps
    Kill all cancer cell neighbours
    If n_gaps>0, cell moves 1 step in a random direction
  Evaluate:
    Mean cancer cells at tmax and t1/2
    Standard deviation at tmax and t1/2
   Mean killing events
   
**Results:**
_
Qualitative observations:_
- Tumour cell automata  follow a logistic growth pattern, with growth rate reaching  plateau at cell plate capacity 
- Introducing immune killing reduces the growth rate (see right), sensitive to agent_ratio
- Mortality rates of c0.5 are required for an equilibrium cell concentration < saturation (@50x50=2500 cells)

Sensitivity analysis (n=3 per instance) reveals:
- A positive correlation between motility and growth constants and cell count at time t
- Increased noise with increased mortality rate
