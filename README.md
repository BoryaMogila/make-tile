## make-tile &mdash; Convert geojson to tile

A highly efficient JavaScript library for **making vector tile from GeoJSON data on the fly**.

Resulting tiles conform to the JSON equivalent
of the [vector tile specification](https://github.com/mapbox/vector-tile-spec/).
To make data rendering and interaction fast, the tiles are simplified,
retaining the minimum level of detail appropriate for each zoom level
(simplifying shapes, filtering out tiny polygons and polylines).


### Usage

```js
// build an tile
const tile = makeTile(geoJSON, options, { x: 149, y: 86, z: 6 });

```

### Options

You can fine-tune the results with an options object,
although the defaults are sensible and work well for most use cases.

```js
const tle = makeTile(
    data,
    {
        maxZoom: 14,  // max zoom to preserve detail on; can't be higher than 24
        tolerance: 3, // simplification tolerance (higher means simpler)
        extent: 4096, // tile extent (both width and height)
        buffer: 64,   // tile buffer on each side
        lineMetrics: false, // whether to enable line metrics tracking for LineString/MultiLineString features
        promoteId: null,    // name of a feature property to promote to feature.id. Cannot be used with `generateId`
        generateId: false,  // whether to generate feature ids. Cannot be used with `promoteId`
        indexMaxPoints: 100000 // max number of points per tile in the index
    },
    {
        x: 149,
        y: 86,
        z: 6
    }
);
```

The `promoteId` and `generateId` options ignore existing `id` values on the feature objects.

### Install

Install using NPM (`npm install make-tile`) or Yarn (`yarn add make-tile`), then:

```js
// import as a ES module
import makeTile from 'make-tile';

// or require in Node / Browserify
const makeTile = require('make-tile');
```
