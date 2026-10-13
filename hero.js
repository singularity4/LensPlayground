/*
 * hero.js
 *
 * The lensed sky at the top of the start page. A field of background galaxies is drawn once in the source plane;
 * every frame, each pixel of the picture is traced back through the lens equation of lensmodel.js to the source
 * plane, so the picture shows the field as lensed by the gold galaxy. The lens follows the pointer.
 * The two small pictures of the start page are drawn in the same way, with a single background galaxy.
 */
(function (global) {
  'use strict';

  const LensModel = global.LensModel;
  const FIELD_HEIGHT = 9;                 // height of the sky shown in the hero, arcsec
  const HERO_LENS = {einsteinRadius: 1.4, e1: 0.06, e2: 0.02, shear1: 0.02, shear2: -0.01};
  const SKY = [6, 10, 20];                // colour of the empty sky (RGB)

  /** Seeded random numbers between 0 and 1 (the "mulberry32" generator), so the field is the same on every visit. */
  function randomGenerator(seed) {
    let state = seed >>> 0;
    return () => {
      state = (state + 0x6D2B79F5) >>> 0;
      let mixed = state;
      mixed = Math.imul(mixed ^ (mixed >>> 15), mixed | 1);
      mixed ^= mixed + Math.imul(mixed ^ (mixed >>> 7), mixed | 61);
      return ((mixed ^ (mixed >>> 14)) >>> 0) / 4294967296;
    };
  }

  /** Background galaxies in the source plane: positions and sizes in arcsec, colour, brightness. */
  function makeGalaxies(width, height, random) {
    const galaxies = [{centreX: 0.12, centreY: 0.05, halfLightRadius: 0.32, sersicIndex: 1.2, axisRatio: 0.6,
                       angle: 0.4, brightness: 1.6, colour: [150, 205, 255]}];   // the galaxy behind the lens at start
    for (let i = 0; i < 110; i++) {
      const reddish = random() < 0.25;
      galaxies.push({
        centreX: (random() - 0.5) * width, centreY: (random() - 0.5) * height,
        halfLightRadius: 0.08 + 0.3 * random() * random(), sersicIndex: 1 + random(),
        axisRatio: 0.35 + 0.65 * random(), angle: Math.PI * random(), brightness: 0.25 + 0.9 * random(),
        colour: reddish ? [255, 214, 170] : [165 + 50 * random(), 205 + 30 * random(), 255],
      });
    }
    return galaxies;
  }

  /** Draws the galaxies into an RGB float image of the source plane (arcsec per pixel = scale). */
  function drawSourcePlane(galaxies, columns, rows, scale) {
    const light = new Float32Array(columns * rows * 3);
    for (const galaxy of galaxies) {
      const [e1, e2] = LensModel.angleAndRatioToEllipticity(galaxy.angle, galaxy.axisRatio);
      const source = {centreX: galaxy.centreX, centreY: galaxy.centreY, halfLightRadius: galaxy.halfLightRadius,
                      sersicIndex: galaxy.sersicIndex, e1, e2};
      const reach = 4 * galaxy.halfLightRadius / Math.sqrt(galaxy.axisRatio);
      const firstColumn = Math.max(0, Math.floor((galaxy.centreX - reach) / scale + columns / 2));
      const lastColumn = Math.min(columns - 1, Math.ceil((galaxy.centreX + reach) / scale + columns / 2));
      const firstRow = Math.max(0, Math.floor((-galaxy.centreY - reach) / scale + rows / 2));
      const lastRow = Math.min(rows - 1, Math.ceil((-galaxy.centreY + reach) / scale + rows / 2));
      for (let row = firstRow; row <= lastRow; row++) {
        for (let column = firstColumn; column <= lastColumn; column++) {
          const x = (column - columns / 2) * scale, y = (rows / 2 - row) * scale;
          const value = galaxy.brightness * LensModel.relativeSersicBrightness(x, y, source);
          const index = 3 * (row * columns + column);
          for (let channel = 0; channel < 3; channel++) light[index + channel] += value * galaxy.colour[channel] / 255;
        }
      }
    }
    return light;
  }

  /** Converts light to display colours with an arcsinh stretch, which keeps faint arcs visible. */
  function toDisplay(light) {
    const pixels = new Uint8ClampedArray(light.length);
    for (let i = 0; i < light.length; i++) {
      const stretched = Math.asinh(6 * light[i]) / Math.asinh(6);
      pixels[i] = Math.min(255, SKY[i % 3] + 255 * stretched);
    }
    return pixels;
  }

  /**
   * Deflection (in pixels, screen orientation) for every offset from the lens centre within +-columns and +-rows,
   * computed once; moving the lens then only shifts where the table is read.
   */
  function deflectionTable(lens, columns, rows, scale) {
    const tableColumns = 2 * columns + 1, tableRows = 2 * rows + 1;
    const alphaColumn = new Float32Array(tableColumns * tableRows);
    const alphaRow = new Float32Array(tableColumns * tableRows);
    for (let row = 0; row < tableRows; row++) {
      for (let column = 0; column < tableColumns; column++) {
        const x = (column - columns) * scale, y = (rows - row) * scale;
        const [alphaX, alphaY] = LensModel.deflection(x, y, lens);
        alphaColumn[row * tableColumns + column] = alphaX / scale;
        alphaRow[row * tableColumns + column] = -alphaY / scale;
      }
    }
    return {alphaColumn, alphaRow, tableColumns};
  }

  /** Lensed picture for a lens centred at pixel (lensColumn, lensRow), by bilinear interpolation in the source plane. */
  function lensImage(output, source, table, columns, rows, lensColumn, lensRow) {
    for (let row = 0; row < rows; row++) {
      for (let column = 0; column < columns; column++) {
        const tableIndex = (row - Math.round(lensRow) + rows) * table.tableColumns + (column - Math.round(lensColumn) + columns);
        const sourceColumn = column - table.alphaColumn[tableIndex];
        const sourceRow = row - table.alphaRow[tableIndex];
        const outputIndex = 4 * (row * columns + column);
        const left = Math.floor(sourceColumn), top = Math.floor(sourceRow);
        if (left < 0 || top < 0 || left >= columns - 1 || top >= rows - 1) {
          output[outputIndex] = SKY[0]; output[outputIndex + 1] = SKY[1]; output[outputIndex + 2] = SKY[2];
        } else {
          const fractionX = sourceColumn - left, fractionY = sourceRow - top;
          const topLeft = 3 * (top * columns + left), topRight = topLeft + 3;
          const bottomLeft = topLeft + 3 * columns, bottomRight = bottomLeft + 3;
          for (let channel = 0; channel < 3; channel++) {
            const upper = (1 - fractionX) * source[topLeft + channel] + fractionX * source[topRight + channel];
            const lower = (1 - fractionX) * source[bottomLeft + channel] + fractionX * source[bottomRight + channel];
            output[outputIndex + channel] = (1 - fractionY) * upper + fractionY * lower;
          }
        }
        output[outputIndex + 3] = 255;
      }
    }
  }

  /** Light of the lens galaxy (foreground, not lensed): a soft gold glow. */
  function drawLensGalaxy(context, column, row, size) {
    const glow = context.createRadialGradient(column, row, 0, column, row, size);
    glow.addColorStop(0, 'rgba(255, 244, 214, 0.95)');
    glow.addColorStop(0.12, 'rgba(240, 214, 160, 0.7)');
    glow.addColorStop(0.45, 'rgba(200, 160, 100, 0.18)');
    glow.addColorStop(1, 'rgba(120, 90, 50, 0)');
    context.globalCompositeOperation = 'screen';
    context.fillStyle = glow;
    context.fillRect(column - size, row - size, 2 * size, 2 * size);
    context.globalCompositeOperation = 'source-over';
  }

  /** Starts the interactive hero in the given canvas. */
  function startHero(canvas) {
    const context = canvas.getContext('2d');
    const reduceMotion = global.matchMedia('(prefers-reduced-motion: reduce)').matches;
    let state = null;

    function setUp() {
      const bounds = canvas.getBoundingClientRect();
      const columns = Math.max(160, Math.round(bounds.width / 2));
      const rows = Math.max(100, Math.round(bounds.height / 2));
      canvas.width = columns;
      canvas.height = rows;
      const scale = FIELD_HEIGHT / rows;
      const galaxies = makeGalaxies(columns * scale, rows * scale, randomGenerator(7));
      // the galaxy behind the lens at start sits under the starting lens position
      const startColumn = 0.66 * columns, startRow = 0.5 * rows;
      galaxies[0].centreX += (startColumn - columns / 2) * scale;
      galaxies[0].centreY += (rows / 2 - startRow) * scale;
      state = {
        columns, rows, scale,
        source: toDisplay(drawSourcePlane(galaxies, columns, rows, scale)),
        table: deflectionTable(HERO_LENS, columns, rows, scale),
        image: context.createImageData(columns, rows),
        lens: [startColumn, startRow], target: [startColumn, startRow],
      };
      draw();
    }

    function draw() {
      const {columns, rows, scale, source, table, image, lens} = state;
      lensImage(image.data, source, table, columns, rows, lens[0], lens[1]);
      context.putImageData(image, 0, 0);
      drawLensGalaxy(context, lens[0], lens[1], 1.6 / scale);
    }

    let animating = false;
    function step() {
      const [targetColumn, targetRow] = state.target;
      const easing = reduceMotion ? 1 : 0.18;
      state.lens[0] += easing * (targetColumn - state.lens[0]);
      state.lens[1] += easing * (targetRow - state.lens[1]);
      draw();
      if (Math.hypot(targetColumn - state.lens[0], targetRow - state.lens[1]) > 0.3) requestAnimationFrame(step);
      else animating = false;
    }

    function moveLensTo(event) {
      const bounds = canvas.getBoundingClientRect();
      state.target = [(event.clientX - bounds.left) / bounds.width * state.columns,
                      (event.clientY - bounds.top) / bounds.height * state.rows];
      if (!animating) { animating = true; requestAnimationFrame(step); }
    }

    canvas.addEventListener('pointermove', moveLensTo);
    canvas.addEventListener('pointerdown', moveLensTo);
    let resizeTimer = null;
    global.addEventListener('resize', () => { clearTimeout(resizeTimer); resizeTimer = setTimeout(setUp, 200); });
    setUp();
  }

  /** Static picture of one lensed galaxy, for the links to the two pages; galaxy.brightness defaults to 1.4. */
  function drawExample(canvas, lens, galaxy) {
    const context = canvas.getContext('2d');
    const columns = canvas.width, rows = canvas.height, scale = 4.2 / rows;
    const [e1, e2] = LensModel.angleAndRatioToEllipticity(galaxy.angle, galaxy.axisRatio);
    const source = {centreX: galaxy.centreX, centreY: galaxy.centreY, halfLightRadius: galaxy.halfLightRadius,
                    sersicIndex: 1.5, e1, e2};
    const image = context.createImageData(columns, rows);
    for (let row = 0; row < rows; row++) {
      for (let column = 0; column < columns; column++) {
        const x = (column - columns / 2) * scale, y = (rows / 2 - row) * scale;
        const [sourceX, sourceY] = LensModel.sourcePosition(x, y, lens);
        const value = Math.asinh(6 * (galaxy.brightness || 1.4) * LensModel.relativeSersicBrightness(sourceX, sourceY, source)) / Math.asinh(6);
        const index = 4 * (row * columns + column);
        image.data[index] = SKY[0] + 150 * value;
        image.data[index + 1] = SKY[1] + 205 * value;
        image.data[index + 2] = SKY[2] + 235 * value;
        image.data[index + 3] = 255;
      }
    }
    context.putImageData(image, 0, 0);
    drawLensGalaxy(context, columns / 2, rows / 2, 0.9 / scale);
  }

  global.Hero = {startHero, drawExample};
})(window);
