import {terser} from 'rollup-plugin-terser';
import buble from '@rollup/plugin-buble';

const config = (file, plugins) => ({
    input: 'src/index.js',
    output: {
        name: 'make-tile',
        format: 'umd',
        indent: false,
        file
    },
    plugins
});

const bubleConfig = {
  transforms: {
    dangerousForOf: true,
  },
  objectAssign: 'Object.assign',
};

export default [
    config('make-tile-dev.js', [buble(bubleConfig)]),
    config('make-tile.js', [terser(), buble(bubleConfig)])
];
