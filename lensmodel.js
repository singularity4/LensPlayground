/*
 * lensmodel.js
 *
 * Lens and light model of LensPlayground. It is the model of the simulated benchmark of EUCLID-SAGE, written in
 * the conventions of lenstronomy 1.14:
 *   lens mass     singular isothermal ellipsoid (lenstronomy "SIE") with external shear ("SHEAR"); with axis
 *                 ratio 1 and no shear it is the singular isothermal sphere (SIS);
 *   source light  elliptical Sersic profile ("SERSIC_ELLIPSE").
 * The Einstein radius follows from the velocity dispersion of the lens galaxy and the redshifts of lens and source
 * in a flat LambdaCDM cosmology (H0 = 70 km/s/Mpc, Omega_m = 0.3). All angles are in arcseconds.
 *
 * A lens is described by {einsteinRadius, e1, e2, shear1, shear2} and a source by
 * {centreX, centreY, halfLightRadius, sersicIndex, e1, e2}, where (e1, e2) are lenstronomy's ellipticity
 * components. tests/compare_with_lenstronomy.py compares every function with lenstronomy and astropy.
 */
(function (global) {
  'use strict';

  const SPEED_OF_LIGHT_KMS = 299792.458;
  const HUBBLE_CONSTANT = 70;              // km/s/Mpc
  const MATTER_DENSITY = 0.3;              // Omega_m of a flat universe
  const ARCSEC_IN_RADIANS = Math.PI / 648000;
  const MEGAPARSEC_IN_METRES = 3.0856775814913673e22;
  const GRAVITATIONAL_CONSTANT = 6.6743e-11;   // m^3 kg^-1 s^-2
  const SPEED_OF_LIGHT_MS = 299792458;
  const SOLAR_MASS_IN_KG = 1.988409870698051e30;
  const SIE_CORE_RADIUS = 1e-10;           // lenstronomy's SIE is its NIE with this core radius

  // ---------------------------------------------------------------------------------------------- cosmology

  /** Comoving distance (Mpc) to redshift z in a flat LambdaCDM universe, by Simpson's rule. */
  function comovingDistance(redshift) {
    const intervals = 400;
    const step = redshift / intervals;
    const inverseHubbleRate = z => 1 / Math.sqrt(MATTER_DENSITY * Math.pow(1 + z, 3) + (1 - MATTER_DENSITY));
    let sum = 0;
    for (let i = 0; i <= intervals; i++) {
      const weight = (i === 0 || i === intervals) ? 1 : (i % 2 === 1 ? 4 : 2);
      sum += weight * inverseHubbleRate(i * step);
    }
    return SPEED_OF_LIGHT_KMS / HUBBLE_CONSTANT * sum * step / 3;
  }

  /** Angular-diameter distances (Mpc) to the lens, to the source, and from the lens to the source. */
  function angularDiameterDistances(lensRedshift, sourceRedshift) {
    const toLens = comovingDistance(lensRedshift);
    const toSource = comovingDistance(sourceRedshift);
    return {
      lens: toLens / (1 + lensRedshift),
      source: toSource / (1 + sourceRedshift),
      lensToSource: (toSource - toLens) / (1 + sourceRedshift),
    };
  }

  /** Einstein radius (arcsec) of an isothermal lens: 4 pi (sigma / c)^2 D_LS / D_S. */
  function einsteinRadius(velocityDispersion, lensRedshift, sourceRedshift) {
    if (sourceRedshift <= lensRedshift) return 0;
    const distances = angularDiameterDistances(lensRedshift, sourceRedshift);
    return 4 * Math.PI * Math.pow(velocityDispersion / SPEED_OF_LIGHT_KMS, 2)
      * distances.lensToSource / distances.source / ARCSEC_IN_RADIANS;
  }

  /** Projected mass (solar masses) inside the Einstein radius: (c^2 / 4G) theta_E^2 D_L D_S / D_LS. */
  function massInsideEinsteinRadius(radius, lensRedshift, sourceRedshift) {
    if (sourceRedshift <= lensRedshift || radius <= 0) return 0;
    const distances = angularDiameterDistances(lensRedshift, sourceRedshift);
    const radiusInRadians = radius * ARCSEC_IN_RADIANS;
    return SPEED_OF_LIGHT_MS * SPEED_OF_LIGHT_MS / (4 * GRAVITATIONAL_CONSTANT)
      * radiusInRadians * radiusInRadians * distances.lens * distances.source / distances.lensToSource
      * MEGAPARSEC_IN_METRES / SOLAR_MASS_IN_KG;
  }

  // -------------------------------------------------------------------------------------------- ellipticity

  /** lenstronomy convention: ellipticity components -> [orientation angle (rad), axis ratio minor/major]. */
  function ellipticityToAngleAndRatio(e1, e2) {
    const modulus = Math.min(Math.hypot(e1, e2), 0.9999);
    return [Math.atan2(e2, e1) / 2, (1 - modulus) / (1 + modulus)];
  }

  /** lenstronomy convention: orientation angle (rad) and axis ratio -> [e1, e2]. */
  function angleAndRatioToEllipticity(angle, axisRatio) {
    const modulus = (1 - axisRatio) / (1 + axisRatio);
    return [modulus * Math.cos(2 * angle), modulus * Math.sin(2 * angle)];
  }

  // --------------------------------------------------------------------------------------------- deflection

  /**
   * Deflection (arcsec) of lenstronomy's SIE, following Keeton & Kochanek (1998). In the frame of the major axis
   * the convergence is kappa = b / (2 sqrt(q^2 x^2 + y^2)); b = theta_E sqrt(q) because lenstronomy's Einstein
   * radius is the geometric mean of the Einstein radii along the major and minor axes.
   */
  function deflectionSIE(x, y, radius, e1, e2) {
    if (radius <= 0) return [0, 0];
    const [angle, axisRatio] = ellipticityToAngleAndRatio(e1, e2);
    const majorAxisRadius = radius / Math.sqrt((1 + axisRatio * axisRatio) / (2 * axisRatio));
    const b = majorAxisRadius * Math.sqrt((1 + axisRatio * axisRatio) / 2);
    const core = SIE_CORE_RADIUS / Math.sqrt(axisRatio);
    const q = axisRatio >= 1 ? 0.99999999 : axisRatio;   // as in lenstronomy, avoids 0/0 for a round lens
    const cosAngle = Math.cos(angle);
    const sinAngle = Math.sin(angle);
    const xMajor = x * cosAngle + y * sinAngle;           // coordinates along the major and minor axes
    const yMinor = -x * sinAngle + y * cosAngle;
    const psi = Math.sqrt(q * q * (core * core + xMajor * xMajor) + yMinor * yMinor);
    const eccentricity = Math.sqrt(1 - q * q);
    const alphaMajor = b / eccentricity * Math.atan(eccentricity * xMajor / (psi + core));
    const alphaMinor = b / eccentricity * Math.atanh(eccentricity * yMinor / (psi + q * q * core));
    return [alphaMajor * cosAngle - alphaMinor * sinAngle, alphaMajor * sinAngle + alphaMinor * cosAngle];
  }

  /** Deflection (arcsec) of an external shear (gamma1, gamma2), lenstronomy "SHEAR". */
  function deflectionShear(x, y, shear1, shear2) {
    return [shear1 * x + shear2 * y, shear2 * x - shear1 * y];
  }

  /** Total deflection (arcsec) of the SIE and the external shear. */
  function deflection(x, y, lens) {
    const isothermal = deflectionSIE(x, y, lens.einsteinRadius, lens.e1 || 0, lens.e2 || 0);
    const shear = deflectionShear(x, y, lens.shear1 || 0, lens.shear2 || 0);
    return [isothermal[0] + shear[0], isothermal[1] + shear[1]];
  }

  /** Lens equation: position in the source plane, beta = theta - alpha(theta). */
  function sourcePosition(x, y, lens) {
    const alpha = deflection(x, y, lens);
    return [x - alpha[0], y - alpha[1]];
  }

  /** Determinant of the lens Jacobian A = d(beta) / d(theta), by central differences. */
  function jacobianDeterminant(x, y, lens) {
    const step = 1e-5;
    const right = deflection(x + step, y, lens), left = deflection(x - step, y, lens);
    const up = deflection(x, y + step, lens), down = deflection(x, y - step, lens);
    const alphaXByX = (right[0] - left[0]) / (2 * step), alphaYByX = (right[1] - left[1]) / (2 * step);
    const alphaXByY = (up[0] - down[0]) / (2 * step), alphaYByY = (up[1] - down[1]) / (2 * step);
    return (1 - alphaXByX) * (1 - alphaYByY) - alphaXByY * alphaYByX;
  }

  // ----------------------------------------------------------------------------- critical curve and caustic

  /**
   * Tangential critical curve, det A = 0, as points [x, y] in the image plane. For an isothermal lens det A < 0 inside
   * the curve and det A > 0 outside, so the curve is found by bisection along rays from the lens centre.
   */
  function criticalCurve(lens, rayCount = 360, maximumRadius = 20) {
    const points = [];
    if (lens.einsteinRadius <= 0) return points;
    for (let i = 0; i < rayCount; i++) {
      const rayAngle = 2 * Math.PI * i / rayCount;
      const cosAngle = Math.cos(rayAngle), sinAngle = Math.sin(rayAngle);
      let inside = 1e-3 * lens.einsteinRadius;
      let outside = maximumRadius;
      if (jacobianDeterminant(inside * cosAngle, inside * sinAngle, lens) >= 0) continue;
      if (jacobianDeterminant(outside * cosAngle, outside * sinAngle, lens) <= 0) continue;
      for (let iteration = 0; iteration < 50; iteration++) {
        const middle = 0.5 * (inside + outside);
        if (jacobianDeterminant(middle * cosAngle, middle * sinAngle, lens) < 0) inside = middle;
        else outside = middle;
      }
      const radius = 0.5 * (inside + outside);
      points.push([radius * cosAngle, radius * sinAngle]);
    }
    return points;
  }

  /** Tangential caustic: the critical curve mapped to the source plane by the lens equation. */
  function caustic(lens, rayCount = 360) {
    return criticalCurve(lens, rayCount).map(([x, y]) => sourcePosition(x, y, lens));
  }

  /**
   * Cut of a singular isothermal lens: the image of the lens centre in the source plane, beta = -alpha(r -> 0).
   * A source outside the cut has no image near the lens centre.
   */
  function cut(lens, rayCount = 360) {
    const points = [];
    if (lens.einsteinRadius <= 0) return points;
    const nearCentre = 1e-7;
    for (let i = 0; i < rayCount; i++) {
      const rayAngle = 2 * Math.PI * i / rayCount;
      const alpha = deflectionSIE(nearCentre * Math.cos(rayAngle), nearCentre * Math.sin(rayAngle),
        lens.einsteinRadius, lens.e1 || 0, lens.e2 || 0);
      points.push([-alpha[0], -alpha[1]]);
    }
    return points;
  }

  // ----------------------------------------------------------------------------------------------- source

  /** lenstronomy's approximation of the Sersic constant b(n). */
  function sersicConstant(sersicIndex) {
    return Math.max(1.9992 * sersicIndex - 0.3271, 1e-5);
  }

  /**
   * Surface brightness of lenstronomy's SERSIC_ELLIPSE with amplitude 1 at the half-light radius, at the
   * position (x, y) in the source plane. The elliptical radius is lenstronomy's: the geometric mean of the distances
   * along the major and minor axes.
   */
  function sersicBrightness(x, y, source) {
    const dx = x - source.centreX, dy = y - source.centreY;
    const e1 = source.e1 || 0, e2 = source.e2 || 0;
    const norm = Math.sqrt(Math.max(Math.abs(1 - e1 * e1 - e2 * e2), 1e-6));
    const xElliptical = ((1 - e1) * dx - e2 * dy) / norm;
    const yElliptical = (-e2 * dx + (1 + e1) * dy) / norm;
    const radius = Math.max(Math.hypot(xElliptical, yElliptical), 1e-4);
    const halfLightRadius = Math.max(source.halfLightRadius, 1e-4);
    const b = sersicConstant(source.sersicIndex);
    return Math.exp(-b * (Math.pow(radius / halfLightRadius, 1 / source.sersicIndex) - 1));
  }

  /** Sersic brightness divided by its central value, between 0 and 1. */
  function relativeSersicBrightness(x, y, source) {
    return sersicBrightness(x, y, source) / Math.exp(sersicConstant(source.sersicIndex));
  }

  const LensModel = {
    comovingDistance, angularDiameterDistances, einsteinRadius, massInsideEinsteinRadius,
    ellipticityToAngleAndRatio, angleAndRatioToEllipticity,
    deflectionSIE, deflectionShear, deflection, sourcePosition, jacobianDeterminant,
    criticalCurve, caustic, cut, sersicConstant, sersicBrightness, relativeSersicBrightness,
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = LensModel;
  else global.LensModel = LensModel;
})(typeof window !== 'undefined' ? window : globalThis);
