# HNNGP tutorial (binary mode choice with rstan)

Tutorial materials for estimating a **hierarchical nearest neighbor Gaussian process (HNNGP)** for binary mode choice using [`rstan`](https://mc-stan.org/users/interfaces/rstan.html). This repository accompanies Villarraga and Daziano (2025).

## Citation

If you use this tutorial, please cite the paper:

> Villarraga, D. F., & Daziano, R. A. (2025). Hierarchical Nearest Neighbor Gaussian Process models for discrete choice: Mode choice in New York City. *Transportation Research Part B: Methodological*, 191, 103132. <https://doi.org/10.1016/j.trb.2024.103132>

A machine-readable citation is in [`CITATION.cff`](CITATION.cff). The code is released under the [MIT License](LICENSE).

Repository: <https://github.com/dfvillarraga2191/hnngp-tutorial>

## Contents

| Path | Description |
|------|-------------|
| `HNNGP_tutorial.Rmd` | Tutorial source (R Markdown) |
| `docs/HNNGP_tutorial.html` | Self-contained HTML (for a website or GitHub Pages) |
| `docs/tutorial.css` | Stylesheet used by the HTML |
| `data/mode_binary.csv` | Example binary mode-choice data |
| `stan/nngp_binary.stan` | Stan model |
| `R/NNMatrix.R` | Nearest-neighbor matrix helpers (`spNNGP`) |
| `R/ess_rhats.R` | ESS (`sns`) and $\widehat{R}$ (`asbio`) diagnostics |
| `scripts/fit_hnngp.R` | End-to-end estimation script |

## Requirements

- R ≥ 4.x with a working C++ toolchain (for Stan)
- Packages: `rstan`, `spNNGP`, `tidyverse`, `HDInterval`, `ggridges`, `loo`, `sns`, `asbio`, `rmarkdown`

```r
install.packages(c(
  "rstan", "spNNGP", "tidyverse", "HDInterval",
  "ggridges", "loo", "sns", "asbio", "rmarkdown"
))
```

## Quick start

From this directory:

```r
source("scripts/fit_hnngp.R")
```

For a faster trial run, uncomment the subsample block inside `scripts/fit_hnngp.R` (or in the Rmd).

## Rebuild the HTML

```r
rmarkdown::render(
  "HNNGP_tutorial.Rmd",
  output_dir = "docs",
  output_file = "HNNGP_tutorial.html"
)
```

The HTML is **self-contained** (`self_contained: true`), so you can upload `docs/HNNGP_tutorial.html` alone, link it from a site, or serve the `docs/` folder with GitHub Pages.

## Website embedding

- **Link / download:** host `docs/HNNGP_tutorial.html` and link to it.
- **iframe:** `<iframe src="path/to/HNNGP_tutorial.html" width="100%" height="900" style="border:0"></iframe>`
- **GitHub Pages:** push this folder and set Pages to `/docs`. The site root opens the tutorial.

## Diagnostics

After fitting, with `fit` in the workspace:

```r
source("R/ess_rhats.R")
```

This prints and returns ESS (summed across chains via `sns::ess`) and Gelman–Rubin $\widehat{R}$ (`asbio::R.hat`) for $\beta$, $\sigma^2$, and $\phi$.

## License / data

Example data are an extract for this tutorial and still contain trip origin and destination coordinates. Household, person, and census-tract identifiers are not included. Adapt paths and priors to your own study as needed.
