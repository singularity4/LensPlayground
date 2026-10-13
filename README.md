# LensPlayground

Interactive galaxy–galaxy strong lensing in the browser. Move a background galaxy behind a lens galaxy to see its
light bend into arcs or a ring.

**Live page:** https://singularity4.github.io/LensPlayground/

| Page | Lens |
|---|---|
| [`circular.html`](circular.html) | singular isothermal sphere (SIS): Einstein ring, arcs, two images |
| [`elliptical.html`](elliptical.html) | singular isothermal ellipsoid (SIE) with external shear: critical curve, caustic, four images |

## Model

LensPlayground uses the same lens and source models as the
[EUCLID-SAGE](https://github.com/singularity4/EUCLID-SAGE) simulations, using lenstronomy 1.14 parameters and conventions:

- **Lens mass:** SIE with external shear (elliptical); the SIS model (circular). 
- **Source light:** elliptical Sérsic profile (`SERSIC_ELLIPSE`), Sérsic index 1.5.
- **Einstein radius:** θ_E = 4π (σ/c)² D_LS / D_S from the velocity dispersion σ of the lens galaxy and the
  redshifts of lens and source, in a flat ΛCDM universe (H0 = 70 km/s/Mpc, Ωm = 0.3). The simulations also show the
  mass inside the Einstein radius.

All angles are in arcseconds. Equations and conventions: [`docs/lens_playground_models.pdf`](docs/lens_playground_models.pdf).

## Files

```
index.html        start page (GitHub Pages)
circular.html     circular lens
elliptical.html   elliptical lens with external shear
lensmodel.js      lens and light model: deflection, critical curve, caustic, Einstein radius, Sersic source
playground.js     drawing and interaction of the two pages
hero.js           example pictures of the start page
playground.css    typeface, colours and layout, light and dark mode
docs/             model note (LaTeX source and PDF)
tests/            comparison of lensmodel.js with lenstronomy
```

## Run locally

Open `index.html` in a browser. No server or installation is needed.

## Check against lenstronomy

```bash
pip install numpy astropy lenstronomy==1.14.2
python tests/compare_with_lenstronomy.py
```

The script runs `lensmodel.js` with Node.js on 20 random configurations and compares deflections, Sérsic brightness,
Einstein radius, mass, critical curve and caustic with lenstronomy and astropy. Deflections agree to 1e-15 arcsec.
