(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? module.exports = factory(require('geojson-vt/src/convert'), require('geojson-vt/src/wrap'), require('geojson-vt/src/transform'), require('geojson-vt/src/tile')) :
typeof define === 'function' && define.amd ? define(['geojson-vt/src/convert', 'geojson-vt/src/wrap', 'geojson-vt/src/transform', 'geojson-vt/src/tile'], factory) :
(global = global || self, global['make-tile'] = factory(global.convert, global.wrap, global.transformTile, global.createTile));
}(this, (function (convert, wrap, transformTile, createTile) { 'use strict';

convert = convert && Object.prototype.hasOwnProperty.call(convert, 'default') ? convert['default'] : convert;
wrap = wrap && Object.prototype.hasOwnProperty.call(wrap, 'default') ? wrap['default'] : wrap;
transformTile = transformTile && Object.prototype.hasOwnProperty.call(transformTile, 'default') ? transformTile['default'] : transformTile;
createTile = createTile && Object.prototype.hasOwnProperty.call(createTile, 'default') ? createTile['default'] : createTile;

var defaultOptions = {
    maxZoom: 14,            // max zoom to preserve detail on
    indexMaxZoom: 5,        // max zoom in the tile index
    indexMaxPoints: 100000, // max number of points per tile in the tile index
    tolerance: 3,           // simplification tolerance (higher means simpler)
    extent: 4096,           // tile extent
    buffer: 64,             // tile buffer on each side
    lineMetrics: false,     // whether to calculate line metrics
    promoteId: null,        // name of a feature property to be promoted to feature.id
    generateId: false,      // whether to generate feature ids. Cannot be used with promoteId
    debug: 0                // logging level (0, 1 or 2)
};

function makeTile(data, tileOptions, pixel) {
    var x = pixel.x;
    var y = pixel.y;
    var z = pixel.z;
    var options = this.options = Object.assign({}, defaultOptions, tileOptions);
    var features = convert(data, options);
    features = wrap(features, options);
    return transformTile(createTile(features, z, x, y, options), options.extent);
}

return makeTile;

})));
