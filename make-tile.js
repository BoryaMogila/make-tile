!function(e,t){"object"==typeof exports&&"undefined"!=typeof module?module.exports=t():"function"==typeof define&&define.amd?define(t):(e=e||self)["make-tile"]=t()}(this,(function(){"use strict";function e(e,t,n,r,i,o){var u=i-n,a=o-r;if(0!==u||0!==a){var l=((e-n)*u+(t-r)*a)/(u*u+a*a);l>1?(n=i,r=o):l>0&&(n+=u*l,r+=a*l)}return(u=e-n)*u+(a=t-r)*a}function t(e,t,r,i){var o={id:void 0===e?null:e,type:t,geometry:r,tags:i,minX:1/0,minY:1/0,maxX:-1/0,maxY:-1/0};return function(e){var t=e.geometry,r=e.type;if("Point"===r||"MultiPoint"===r||"LineString"===r)n(e,t);else if("Polygon"===r||"MultiLineString"===r)for(var i=0;i<t.length;i++)n(e,t[i]);else if("MultiPolygon"===r)for(i=0;i<t.length;i++)for(var o=0;o<t[i].length;o++)n(e,t[i][o])}(o),o}function n(e,t){for(var n=0;n<t.length;n+=3)e.minX=Math.min(e.minX,t[n]),e.minY=Math.min(e.minY,t[n+1]),e.maxX=Math.max(e.maxX,t[n]),e.maxY=Math.max(e.maxY,t[n+1])}function r(e,n,a,l){if(n.geometry){var s=n.geometry.coordinates,f=n.geometry.type,g=Math.pow(a.tolerance/((1<<a.maxZoom)*a.extent),2),h=[],m=n.id;if(a.promoteId?m=n.properties[a.promoteId]:a.generateId&&(m=l||0),"Point"===f)i(s,h);else if("MultiPoint"===f)for(var p=0;p<s.length;p++)i(s[p],h);else if("LineString"===f)o(s,h,g,!1);else if("MultiLineString"===f){if(a.lineMetrics){for(p=0;p<s.length;p++)h=[],o(s[p],h,g,!1),e.push(t(m,"LineString",h,n.properties));return}u(s,h,g,!1)}else if("Polygon"===f)u(s,h,g,!0);else{if("MultiPolygon"!==f){if("GeometryCollection"===f){for(p=0;p<n.geometry.geometries.length;p++)r(e,{id:m,geometry:n.geometry.geometries[p],properties:n.properties},a,l);return}throw new Error("Input data is not a valid GeoJSON object.")}for(p=0;p<s.length;p++){var v=[];u(s[p],v,g,!0),h.push(v)}}e.push(t(m,f,h,n.properties))}}function i(e,t){t.push(a(e[0])),t.push(l(e[1])),t.push(0)}function o(t,n,r,i){for(var o,u,s=0,f=0;f<t.length;f++){var g=a(t[f][0]),h=l(t[f][1]);n.push(g),n.push(h),n.push(0),f>0&&(s+=i?(o*h-g*u)/2:Math.sqrt(Math.pow(g-o,2)+Math.pow(h-u,2))),o=g,u=h}var m=n.length-3;n[2]=1,function t(n,r,i,o){for(var u,a=o,l=i-r>>1,s=i-r,f=n[r],g=n[r+1],h=n[i],m=n[i+1],p=r+3;p<i;p+=3){var v=e(n[p],n[p+1],f,g,h,m);if(v>a)u=p,a=v;else if(v===a){var c=Math.abs(p-l);c<s&&(u=p,s=c)}}a>o&&(u-r>3&&t(n,r,u,o),n[u+2]=a,i-u>3&&t(n,u,i,o))}(n,0,m,r),n[m+2]=1,n.size=Math.abs(s),n.start=0,n.end=n.size}function u(e,t,n,r){for(var i=0;i<e.length;i++){var u=[];o(e[i],u,n,r),t.push(u)}}function a(e){return e/360+.5}function l(e){var t=Math.sin(e*Math.PI/180),n=.5-.25*Math.log((1+t)/(1-t))/Math.PI;return n<0?0:n>1?1:n}function s(e,n,r,i,o,u,a,l){if(i/=n,u>=(r/=n)&&a<i)return e;if(a<r||u>=i)return null;for(var s=[],h=0;h<e.length;h++){var p=e[h],v=p.geometry,c=p.type,y=0===o?p.minX:p.minY,M=0===o?p.maxX:p.maxY;if(y>=r&&M<i)s.push(p);else if(!(M<r||y>=i)){var d=[];if("Point"===c||"MultiPoint"===c)f(v,d,r,i,o);else if("LineString"===c)g(v,d,r,i,o,!1,l.lineMetrics);else if("MultiLineString"===c)m(v,d,r,i,o,!1);else if("Polygon"===c)m(v,d,r,i,o,!0);else if("MultiPolygon"===c)for(var x=0;x<v.length;x++){var P=[];m(v[x],P,r,i,o,!0),P.length&&d.push(P)}if(d.length){if(l.lineMetrics&&"LineString"===c){for(x=0;x<d.length;x++)s.push(t(p.id,c,d[x],p.tags));continue}"LineString"!==c&&"MultiLineString"!==c||(1===d.length?(c="LineString",d=d[0]):c="MultiLineString"),"Point"!==c&&"MultiPoint"!==c||(c=3===d.length?"Point":"MultiPoint"),s.push(t(p.id,c,d,p.tags))}}}return s.length?s:null}function f(e,t,n,r,i){for(var o=0;o<e.length;o+=3){var u=e[o+i];u>=n&&u<=r&&(t.push(e[o]),t.push(e[o+1]),t.push(e[o+2]))}}function g(e,t,n,r,i,o,u){for(var a,l,s=h(e),f=0===i?v:c,g=e.start,m=0;m<e.length-3;m+=3){var y=e[m],M=e[m+1],d=e[m+2],x=e[m+3],P=e[m+4],S=0===i?y:M,L=0===i?x:P,X=!1;u&&(a=Math.sqrt(Math.pow(y-x,2)+Math.pow(M-P,2))),S<n?L>n&&(l=f(s,y,M,x,P,n),u&&(s.start=g+a*l)):S>r?L<r&&(l=f(s,y,M,x,P,r),u&&(s.start=g+a*l)):p(s,y,M,d),L<n&&S>=n&&(l=f(s,y,M,x,P,n),X=!0),L>r&&S<=r&&(l=f(s,y,M,x,P,r),X=!0),!o&&X&&(u&&(s.end=g+a*l),t.push(s),s=h(e)),u&&(g+=a)}var Y=e.length-3;y=e[Y],M=e[Y+1],d=e[Y+2],(S=0===i?y:M)>=n&&S<=r&&p(s,y,M,d),Y=s.length-3,o&&Y>=3&&(s[Y]!==s[0]||s[Y+1]!==s[1])&&p(s,s[0],s[1],s[2]),s.length&&t.push(s)}function h(e){var t=[];return t.size=e.size,t.start=e.start,t.end=e.end,t}function m(e,t,n,r,i,o){for(var u=0;u<e.length;u++)g(e[u],t,n,r,i,o,!1)}function p(e,t,n,r){e.push(t),e.push(n),e.push(r)}function v(e,t,n,r,i,o){var u=(o-t)/(r-t);return e.push(o),e.push(n+(i-n)*u),e.push(1),u}function c(e,t,n,r,i,o){var u=(o-n)/(i-n);return e.push(t+(r-t)*u),e.push(o),e.push(1),u}function y(e,n){for(var r=[],i=0;i<e.length;i++){var o,u=e[i],a=u.type;if("Point"===a||"MultiPoint"===a||"LineString"===a)o=M(u.geometry,n);else if("MultiLineString"===a||"Polygon"===a){o=[];for(var l=0;l<u.geometry.length;l++)o.push(M(u.geometry[l],n))}else if("MultiPolygon"===a)for(o=[],l=0;l<u.geometry.length;l++){for(var s=[],f=0;f<u.geometry[l].length;f++)s.push(M(u.geometry[l][f],n));o.push(s)}r.push(t(u.id,a,o,u.tags))}return r}function M(e,t){var n=[];n.size=e.size,void 0!==e.start&&(n.start=e.start,n.end=e.end);for(var r=0;r<e.length;r+=3)n.push(e[r]+t,e[r+1],e[r+2]);return n}function d(e,t,n,r,i,o){return[Math.round(n*(e*r-i)),Math.round(n*(t*r-o))]}function x(e,t,n,r){var i=t.geometry,o=t.type,u=[];if("Point"===o||"MultiPoint"===o)for(var a=0;a<i.length;a+=3)u.push(i[a]),u.push(i[a+1]),e.numPoints++,e.numSimplified++;else if("LineString"===o)P(u,i,e,n,!1,!1);else if("MultiLineString"===o||"Polygon"===o)for(a=0;a<i.length;a++)P(u,i[a],e,n,"Polygon"===o,0===a);else if("MultiPolygon"===o)for(var l=0;l<i.length;l++){var s=i[l];for(a=0;a<s.length;a++)P(u,s[a],e,n,!0,0===a)}if(u.length){var f=t.tags||null;if("LineString"===o&&r.lineMetrics){for(var g in f={},t.tags)f[g]=t.tags[g];f.mapbox_clip_start=i.start/i.size,f.mapbox_clip_end=i.end/i.size}var h={geometry:u,type:"Polygon"===o||"MultiPolygon"===o?3:"LineString"===o||"MultiLineString"===o?2:1,tags:f};null!==t.id&&(h.id=t.id),e.features.push(h)}}function P(e,t,n,r,i,o){var u=r*r;if(r>0&&t.size<(i?u:r))n.numPoints+=t.length/3;else{for(var a=[],l=0;l<t.length;l+=3)(0===r||t[l+2]>u)&&(n.numSimplified++,a.push(t[l]),a.push(t[l+1])),n.numPoints++;i&&function(e,t){for(var n=0,r=0,i=e.length,o=i-2;r<i;o=r,r+=2)n+=(e[r]-e[o])*(e[r+1]+e[o+1]);if(n>0===t)for(r=0,i=e.length;r<i/2;r+=2){var u=e[r],a=e[r+1];e[r]=e[i-2-r],e[r+1]=e[i-1-r],e[i-2-r]=u,e[i-1-r]=a}}(a,o),e.push(a)}}var S={maxZoom:14,indexMaxZoom:5,indexMaxPoints:1e5,tolerance:3,extent:4096,buffer:64,lineMetrics:!1,promoteId:null,generateId:!1,debug:0};return function(e,t,n){var i=n.x,o=n.y,u=n.z,a=Object.assign({},S,t),l=function(e,t){var n=[];if("FeatureCollection"===e.type)for(var i=0;i<e.features.length;i++)r(n,e.features[i],t,i);else"Feature"===e.type?r(n,e,t):r(n,{geometry:e},t);return n}(e,a);return function(e,t){if(e.transformed)return e;var n,r,i,o=1<<e.z,u=e.x,a=e.y;for(n=0;n<e.features.length;n++){var l=e.features[n],s=l.geometry,f=l.type;if(l.geometry=[],1===f)for(r=0;r<s.length;r+=2)l.geometry.push(d(s[r],s[r+1],t,o,u,a));else for(r=0;r<s.length;r++){var g=[];for(i=0;i<s[r].length;i+=2)g.push(d(s[r][i],s[r][i+1],t,o,u,a));l.geometry.push(g)}}return e.transformed=!0,e}(function(e,t,n,r,i){for(var o=t===i.maxZoom?0:i.tolerance/((1<<t)*i.extent),u={features:[],numPoints:0,numSimplified:0,numFeatures:0,source:null,x:n,y:r,z:t,transformed:!1,minX:2,minY:1,maxX:-1,maxY:0},a=0;a<e.length;a++){u.numFeatures++,x(u,e[a],o,i);var l=e[a].minX,s=e[a].minY,f=e[a].maxX,g=e[a].maxY;l<u.minX&&(u.minX=l),s<u.minY&&(u.minY=s),f>u.maxX&&(u.maxX=f),g>u.maxY&&(u.maxY=g)}return u}(l=function(e,t){var n=t.buffer/t.extent,r=e,i=s(e,1,-1-n,n,0,-1,2,t),o=s(e,1,1-n,2+n,0,-1,2,t);return(i||o)&&(r=s(e,1,-n,1+n,0,-1,2,t)||[],i&&(r=y(i,1).concat(r)),o&&(r=r.concat(y(o,-1)))),r}(l,a),u,i,o,a),a.extent)}}));
