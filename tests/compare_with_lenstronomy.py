"""Compare lensmodel.js with lenstronomy and astropy.

The script evaluates lensmodel.js with Node.js for 20 random configurations of an SIE with external shear and an
elliptical Sersic source, and compares the result with lenstronomy and astropy: deflection, determinant of the lens
Jacobian, Sersic brightness, Einstein radius from the velocity dispersion, mass inside the Einstein radius, and the
critical curve and caustic.

Requirements: Node.js, numpy, astropy and lenstronomy (tested with lenstronomy 1.14.2).
Usage: python tests/compare_with_lenstronomy.py
"""
import json
import pathlib
import subprocess
import tempfile

import astropy.constants as const
import astropy.units as u
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from lenstronomy.Cosmo.lens_cosmo import LensCosmo
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LightModel.light_model import LightModel

LENS_MODEL_FILE = pathlib.Path(__file__).resolve().parents[1] / "lensmodel.js"
NODE_PROGRAM = """
const LensModel = require(process.argv[2]);
const fs = require('fs');
const configurations = JSON.parse(fs.readFileSync(process.argv[3]));
const results = configurations.map(c => ({
  deflection: c.points.map(([x, y]) => LensModel.deflection(x, y, c.lens)),
  jacobianDeterminant: c.points.map(([x, y]) => LensModel.jacobianDeterminant(x, y, c.lens)),
  sersicBrightness: c.points.map(([x, y]) => LensModel.sersicBrightness(x, y, c.source)),
  einsteinRadius: LensModel.einsteinRadius(c.velocityDispersion, c.lensRedshift, c.sourceRedshift),
  massInsideOneArcsec: LensModel.massInsideEinsteinRadius(1.0, c.lensRedshift, c.sourceRedshift),
  criticalCurve: LensModel.criticalCurve(c.lens, 720),
  caustic: LensModel.caustic(c.lens, 720),
}));
fs.writeFileSync(process.argv[4], JSON.stringify(results));
"""


def random_configurations(count=20, seed=3):
    """Random lenses, sources, redshifts and positions in the image plane.

    The first lens is a singular isothermal sphere.
    """
    rng = np.random.default_rng(seed)
    configurations = []
    for index in range(count):
        lens = {"einsteinRadius": rng.uniform(0.3, 2.5), "e1": rng.normal(0, 0.15), "e2": rng.normal(0, 0.15),
                "shear1": rng.normal(0, 0.05), "shear2": rng.normal(0, 0.05)}
        if index == 0:
            lens.update(e1=0.0, e2=0.0, shear1=0.0, shear2=0.0)
        source = {"centreX": rng.uniform(-0.5, 0.5), "centreY": rng.uniform(-0.5, 0.5),
                  "halfLightRadius": rng.uniform(0.05, 0.5), "sersicIndex": rng.uniform(1, 2.5),
                  "e1": rng.normal(0, 0.2), "e2": rng.normal(0, 0.2)}
        configurations.append({
            "lens": {key: float(value) for key, value in lens.items()},
            "source": {key: float(value) for key, value in source.items()},
            "points": rng.uniform(-3, 3, size=(50, 2)).tolist(),
            "velocityDispersion": float(rng.uniform(150, 300)),
            "lensRedshift": float(rng.uniform(0.1, 1.0)),
            "sourceRedshift": float(rng.uniform(1.1, 4.0)),
        })
    return configurations


def evaluate_javascript(configurations):
    """Runs lensmodel.js with Node.js and returns its results for every configuration."""
    with tempfile.TemporaryDirectory() as folder:
        folder = pathlib.Path(folder)
        (folder / "program.js").write_text(NODE_PROGRAM)
        (folder / "configurations.json").write_text(json.dumps(configurations))
        subprocess.run(["node", str(folder / "program.js"), str(LENS_MODEL_FILE), str(folder / "configurations.json"),
                        str(folder / "results.json")], check=True)
        return json.loads((folder / "results.json").read_text())


def jacobian_determinant(lens_model, x, y, lens_keywords):
    """Determinant of the lens Jacobian from lenstronomy's Hessian."""
    psi_xx, psi_xy, psi_yx, psi_yy = lens_model.hessian(x, y, lens_keywords)
    return (1 - psi_xx) * (1 - psi_yy) - psi_xy * psi_yx


def main():
    configurations = random_configurations()
    javascript = evaluate_javascript(configurations)
    lens_model = LensModel(["SIE", "SHEAR"])
    light_model = LightModel(["SERSIC_ELLIPSE"])
    cosmology = FlatLambdaCDM(H0=70, Om0=0.3)
    differences = {name: 0.0 for name in ["deflection (arcsec)", "det A", "Sersic brightness (relative)",
                                          "Einstein radius (relative)", "mass (relative)",
                                          "det A on the critical curve", "caustic (arcsec)"]}

    def record(name, value):
        differences[name] = max(differences[name], float(value))

    for configuration, result in zip(configurations, javascript, strict=True):
        lens, source = configuration["lens"], configuration["source"]
        lens_keywords = [{"theta_E": lens["einsteinRadius"], "e1": lens["e1"], "e2": lens["e2"],
                          "center_x": 0, "center_y": 0},
                         {"gamma1": lens["shear1"], "gamma2": lens["shear2"]}]
        source_keywords = [{"amp": 1, "R_sersic": source["halfLightRadius"], "n_sersic": source["sersicIndex"],
                            "e1": source["e1"], "e2": source["e2"],
                            "center_x": source["centreX"], "center_y": source["centreY"]}]
        x, y = np.array(configuration["points"]).T

        deflection = np.c_[lens_model.alpha(x, y, lens_keywords)]
        record("deflection (arcsec)", np.abs(np.array(result["deflection"]) - deflection).max())
        record("det A", np.abs(np.array(result["jacobianDeterminant"])
                               - jacobian_determinant(lens_model, x, y, lens_keywords)).max())
        brightness = light_model.surface_brightness(x, y, source_keywords)
        record("Sersic brightness (relative)", (np.abs(np.array(result["sersicBrightness"]) - brightness)
                                                / brightness).max())

        lens_redshift, source_redshift = configuration["lensRedshift"], configuration["sourceRedshift"]
        lens_cosmology = LensCosmo(z_lens=lens_redshift, z_source=source_redshift, cosmo=cosmology)
        einstein_radius = lens_cosmology.sis_sigma_v2theta_E(configuration["velocityDispersion"])
        record("Einstein radius (relative)", abs(result["einsteinRadius"] - einstein_radius) / einstein_radius)
        distance_lens = cosmology.angular_diameter_distance(lens_redshift)
        distance_source = cosmology.angular_diameter_distance(source_redshift)
        distance_lens_source = cosmology.angular_diameter_distance_z1z2(lens_redshift, source_redshift)
        mass = (const.c ** 2 / (4 * const.G) * (1 * u.arcsec).to(u.rad).value ** 2
                * distance_lens * distance_source / distance_lens_source).to(u.Msun).value
        record("mass (relative)", abs(result["massInsideOneArcsec"] - mass) / mass)

        critical_x, critical_y = np.array(result["criticalCurve"]).T
        record("det A on the critical curve",
               np.abs(jacobian_determinant(lens_model, critical_x, critical_y, lens_keywords)).max())
        caustic = np.c_[lens_model.ray_shooting(critical_x, critical_y, lens_keywords)]
        record("caustic (arcsec)", np.abs(np.array(result["caustic"]) - caustic).max())

    for name, value in differences.items():
        print(f"largest difference, {name:30s} {value:.1e}")
    assert differences["deflection (arcsec)"] < 1e-10
    assert differences["Sersic brightness (relative)"] < 1e-10
    assert differences["Einstein radius (relative)"] < 1e-6 and differences["mass (relative)"] < 1e-6
    assert differences["det A on the critical curve"] < 1e-3 and differences["caustic (arcsec)"] < 1e-10
    print("lensmodel.js agrees with lenstronomy and astropy.")


if __name__ == "__main__":
    main()
