import convert from './convert';
import wrap from './wrap';
import transformTile from './transform';
import createTile from './tile';

const defaultOptions = {
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

export default function makeTile(data, tileOptions, pixel) {
    const { x, y, z } = pixel;
    const options = Object.assign({}, defaultOptions, tileOptions);
    let features = convert(data, options);
    features = wrap(features, options);
    return transformTile(createTile(features, z, x, y, options), options.extent);
}
