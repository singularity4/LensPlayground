/*
 * playground.js
 *
 * Drawing and interaction shared by the two pages of LensPlayground. The large panel is the image plane, the
 * lensed source as a telescope would see it; the inset is the source plane, the position of the source without
 * the lens. Both use the model of lensmodel.js. Angles are in arcseconds; the y axis points up, as on the sky.
 */
(function (global) {
  'use strict';

  const LensModel = global.LensModel;   // from lensmodel.js, loaded before this file
  const SVG_NAMESPACE = 'http://www.w3.org/2000/svg';
  const IMAGE_SCALE = 70;                 // image plane: SVG units per arcsecond
  const IMAGE_CENTRE = [310, 250];        // lens centre in the image plane
  const INSET_SCALE = 40;                 // inset of the source plane: SVG units per arcsecond
  const INSET_CORNER = [458, 18];         // top-left corner of the inset in the SVG of the image plane
  const INSET_CENTRE = [72, 78];          // optical axis in the coordinates of the inset
  const SAMPLE_COUNT = 30000;             // points in the image plane at which the lensed source is evaluated
  const SAMPLE_RADIUS = 4.4;              // radius (arcsec) of the sampled disc
  const GOLDEN_ANGLE = Math.PI * (3 - Math.sqrt(5));   // spreads the samples evenly over the disc
  const FAINTEST_DRAWN = 0.04;            // samples fainter than 4% of the central source brightness are not drawn
  const BRIGHT_LEVEL = 0.35;              // samples brighter than 35% are drawn in the lighter colour
  const SOURCE_SERSIC_INDEX = 1.5;        // fixed Sersic index of the source (the simulations draw 1 to 2.5)

  /** Faint background stars, for decoration only. */
  function drawStars(group) {
    for (let i = 0; i < 80; i++) {
      const star = document.createElementNS(SVG_NAMESPACE, 'circle');
      star.setAttribute('cx', (i * 83) % 610 + 5);
      star.setAttribute('cy', (i * 137) % 475 + 5);
      star.setAttribute('r', i % 9 === 0 ? 1.1 : 0.45);
      star.setAttribute('fill', '#d7dce5');
      star.setAttribute('opacity', 0.2 + (i % 5) * 0.1);
      group.appendChild(star);
    }
  }

  /**
   * Draws the lensed source: each sampled point in the image plane is mapped to the source plane by the lens equation
   * and drawn with the source brightness there (lensing conserves surface brightness).
   */
  function drawLensedSource(group, lens, source) {
    group.replaceChildren();
    for (let i = 0; i < SAMPLE_COUNT; i++) {
      const sampleAngle = i * GOLDEN_ANGLE;
      const sampleRadius = Math.sqrt((i + 0.5) / SAMPLE_COUNT) * SAMPLE_RADIUS;
      const x = sampleRadius * Math.cos(sampleAngle);
      const y = sampleRadius * Math.sin(sampleAngle);
      const [sourceX, sourceY] = LensModel.sourcePosition(x, y, lens);
      const brightness = LensModel.relativeSersicBrightness(sourceX, sourceY, source);
      if (brightness <= FAINTEST_DRAWN) continue;
      const point = document.createElementNS(SVG_NAMESPACE, 'circle');
      point.setAttribute('cx', IMAGE_CENTRE[0] + x * IMAGE_SCALE);
      point.setAttribute('cy', IMAGE_CENTRE[1] - y * IMAGE_SCALE);
      point.setAttribute('r', 1.45);
      point.setAttribute('fill', brightness > BRIGHT_LEVEL ? '#d8f1ff' : '#78bfe0');
      point.setAttribute('opacity', Math.min(0.9, 0.1 + Math.sqrt(brightness)));
      group.appendChild(point);
    }
  }

  /** SVG path through the points [x, y] (arcsec), closed. */
  function closedPath(points, centre, scale) {
    if (points.length === 0) return '';
    const corners = points.map(([x, y]) => (centre[0] + x * scale).toFixed(2) + ',' + (centre[1] - y * scale).toFixed(2));
    return 'M ' + corners.join(' L ') + ' Z';
  }

  /** Path of a curve in the image plane. */
  function imagePlanePath(points) {
    return closedPath(points, IMAGE_CENTRE, IMAGE_SCALE);
  }

  /** Path of a curve in the inset of the source plane. */
  function sourcePlanePath(points) {
    return closedPath(points, INSET_CENTRE, INSET_SCALE);
  }

  /** Radius in the image plane, in SVG units. */
  function imagePlaneLength(arcseconds) {
    return arcseconds * IMAGE_SCALE;
  }

  /** Radius in the inset of the source plane, in SVG units. */
  function sourcePlaneLength(arcseconds) {
    return arcseconds * INSET_SCALE;
  }

  /** Places the source marker of the inset: an ellipse with the half-light radius and axis ratio of the source. */
  function placeSourceMarker(marker, source, axisRatio) {
    marker.setAttribute('cx', INSET_CENTRE[0] + source.centreX * INSET_SCALE);
    marker.setAttribute('cy', INSET_CENTRE[1] - source.centreY * INSET_SCALE);
    marker.setAttribute('rx', Math.max(2.5, source.halfLightRadius / Math.sqrt(axisRatio) * INSET_SCALE));
    marker.setAttribute('ry', Math.max(2, source.halfLightRadius * Math.sqrt(axisRatio) * INSET_SCALE));
  }

  /** Lets the reader drag the source marker; the sliders of the source position follow the marker. */
  function enableSourceDragging(svg, marker, sliderX, sliderY, onMove) {
    const limitX = +sliderX.max, limitY = +sliderY.max;
    const axisX = INSET_CORNER[0] + INSET_CENTRE[0], axisY = INSET_CORNER[1] + INSET_CENTRE[1];
    let dragging = false;
    marker.addEventListener('pointerdown', event => {
      dragging = true;
      marker.setPointerCapture(event.pointerId);
    });
    marker.addEventListener('pointermove', event => {
      if (!dragging) return;
      const pointer = svg.createSVGPoint();
      pointer.x = event.clientX;
      pointer.y = event.clientY;
      const position = pointer.matrixTransform(svg.getScreenCTM().inverse());
      sliderX.value = Math.max(-limitX, Math.min(limitX, (position.x - axisX) / INSET_SCALE));
      sliderY.value = Math.max(-limitY, Math.min(limitY, (axisY - position.y) / INSET_SCALE));
      onMove();
    });
    marker.addEventListener('pointerup', () => { dragging = false; });
    marker.addEventListener('pointercancel', () => { dragging = false; });
  }

  /** Keeps the source at least 0.05 in redshift behind the lens. */
  function keepSourceBehindLens(lensRedshiftSlider, sourceRedshiftSlider) {
    const minimum = +lensRedshiftSlider.value + 0.05;
    if (+sourceRedshiftSlider.value < minimum) sourceRedshiftSlider.value = Math.min(+sourceRedshiftSlider.max, minimum);
  }

  /** Angle with sign, e.g. "+0.08″". */
  function formatAngle(arcseconds) {
    return (arcseconds > 0 ? '+' : '') + arcseconds.toFixed(2) + '″';
  }

  /** Mass in solar masses in scientific notation, e.g. "3.2 × 10¹¹ M☉". */
  function formatMass(solarMasses) {
    if (!(solarMasses > 0)) return '–';
    const exponent = Math.floor(Math.log10(solarMasses));
    const superscript = String(exponent).replace(/[0-9-]/g, digit => '⁰¹²³⁴⁵⁶⁷⁸⁹⁻'['0123456789-'.indexOf(digit)]);
    return (solarMasses / Math.pow(10, exponent)).toFixed(1) + ' × 10' + superscript + ' M☉';
  }

  global.Playground = {
    SOURCE_SERSIC_INDEX, drawStars, drawLensedSource, imagePlanePath, sourcePlanePath, imagePlaneLength,
    sourcePlaneLength, placeSourceMarker, enableSourceDragging, keepSourceBehindLens, formatAngle, formatMass,
  };
})(window);
