# Does a biased fossil record fool ancestral state estimation?

A simulation study testing whether **state-dependent fossil preservation bias**
distorts the accuracy of **ancestral state estimation (ASE)** for a binary trait.

## Why this matters

When biologists want to know what an extinct ancestor was probably like — did
it have fur or scales, was it a burrower or a swimmer — they usually infer it
statistically, using the traits of living descendants and a model of how traits
evolve on a family tree. This is called **ancestral state estimation**.

Adding fossils to that family tree is known to make these estimates more
accurate, because fossils give a direct, older glimpse of the trait rather than
requiring inference all the way back from the present day.

But there's a catch: **fossils are not a random sample of the past.** Some
traits or lifestyles fossilise far more easily than others — an organism living
somewhere with fine sediment and no scavengers is much more likely to leave a
fossil than one living somewhere that destroys remains quickly. If a trait is
correlated with how *preservable* an organism is, the fossil record itself
becomes a biased witness. A biologist using biased fossils might not realise
their "improved" estimate is actually being pulled in a systematic direction.

**This project asks: how much does that bias actually matter?** Specifically,
does *heterogeneous* preservation (one trait state preserving very differently
to the other) make ancestral state estimates measurably worse than
*homogeneous* preservation (both states preserving at similar rates), even
when the total amount of fossil evidence is held constant?

This is a simulation study rather than a study of real fossils, because with
simulations we always know the ground truth — the *actual* history of the
trait — and can therefore directly measure how wrong an estimation method is.
That's not possible with real data, where the true ancestral states are
exactly what we're trying to find out.

## What the simulation actually does

For each simulated run, the pipeline:

1. **Simulates a family tree** (birth–death process, ~50+ living species)
2. **Evolves a binary trait** ("A" or "B") along that tree under a simple
   random-switching model
3. **Generates a fossil record**, where fossils of a lineage currently in
   state A are preserved at rate `psi_A`, and lineages in state B at rate
   `psi_B` — these two rates can be set equal (no bias) or very different
   (strong bias)
4. **Builds the tree a palaeontologist would actually see**: the living
   species plus only the fossils that happened to be preserved
5. **Estimates the ancestral states** from that fossil-inclusive tree, using
   a standard statistical model (equal-rates Markov model)
6. **Scores the estimate** against the true, simulated history — something
   only possible because this is a simulation, not real data

This is repeated 500 times for each of 16 preservation scenarios (5 levels of
preservation *heterogeneity* × 3 levels of overall preservation *intensity*,
plus a no-fossil control) — 8,000 simulated runs in total.

## Key findings

- Preservation heterogeneity (state A and B fossilising at very different
  rates) makes ancestral state estimates **measurably and consistently
  worse**, and the effect gets steadily larger as the bias increases
- **The total amount of fossil evidence matters more than the bias does** —
  overall preservation intensity explains roughly 4–5× more of the variation
  in accuracy than heterogeneity does
- Even biased fossil evidence is still much better than none: every
  fossil scenario tested clearly outperformed the no-fossil baseline

## Repository structure

```
ase-fossil-simulation/
├── README.md
├── .gitignore
└── R/
    ├── 01_simulate_fossils.R    # runs the 8,000 simulations, saves a CSV
    └── 02_analyze_results.R     # reads that CSV, makes tables + figures
```

The scripts are numbered because they're meant to run in that order —
`01` produces the raw data, `02` analyses it. Everything is deliberately
written as **two separate steps** so you can re-explore the results (`02`)
without needing to re-run hours of simulation (`01`) every time.

## Requirements

R (4.x) with the following packages:

```r
install.packages(c(
  "phytools", "phangorn", "TreeSim", "FossilSim", "paleotree",
  "dplyr", "tidyr", "ggplot2", "emmeans"
))
```

## How to use this

- **To just look at the results**: you'll need a copy of the simulation
  output CSV (not included in this repo — see note below) and can go
  straight to `02_analyze_results.R`.
- **To reproduce everything from scratch**: run `01_simulate_fossils.R`
  first. Be warned — simulating and analysing 8,000 individual phylogenies
  takes a long time (realistically several hours). It saves progress to CSV
  every 100 successful runs, so it's safe to interrupt and resume.

> **Note on data files:** the raw simulation results (CSV, several MB) and
> the figures/tables `02_analyze_results.R` produces aren't committed to this
> repository, to keep it lightweight. Run `01_simulate_fossils.R` to
> regenerate them, or add your own copy of the results CSV to the project
> folder before running `02_analyze_results.R`.

## Status

This is an active MSc dissertation project — the analysis and write-up are
still being revised based on supervisor feedback.
