/**
* Espruino Ublox NEO6MV2 GPS + PCD8544 LCD DRIVER (NOKIA 5110)
* v0.1
* Convert between GPS WGS84 to OSGB36 OS grid reference
*/

var nsum = 0, esum = 0, navg = 'searching', eavg = 'searching', g;
var nelmt = [];
var eelmt = [];

/*  
* Ordnance Survey Grid Reference functions  (c) Chris Veness 2005-2014
*   - www.movable-type.co.uk/scripts/gridref.js
*   - www.movable-type.co.uk/scripts/latlon-gridref.html
*/  
/**
 * Creates a OsGridRef object
 *
 * @constructor
 * @classdesc Convert OS grid references to/from OSGB latitude/longitude points
 * @requires LatLonE, GeoParams (from latlon-ellipse.js)
 *
 * @param {Number} easting - Easting in metres from OS false origin
 * @param {Number} northing - Northing in metres from OS false origin
 */
function OsGridRef(easting, northing) {
    this.easting = Math.floor(Number(easting));
    this.northing = Math.floor(Number(northing));
}


/**
 * Convert (OSGB36) latitude/longitude to Ordnance Survey grid reference easting/northing coordinate
 *
 * @param   {LatLonE} point - OSGB36 latitude/longitude
 * @returns {OsGridRef} OS Grid Reference easting/northing
 * @throws  Error if datum of point is not OSGB36
 */
OsGridRef.latLongToOsGrid = function(point) {
    //if (point.datum != GeoParams.datum.OSGB36) throw new Error('Can only convert OSGB36 point to OsGrid');
	var wgs84 = new LatLonE(point.lat, point.lon, GeoParams.datum.WGS84);
    var osgb = wgs84.convertDatum(GeoParams.datum.OSGB36);
    var varphi = osgb.lat.toRadians();
    var lamda = osgb.lon.toRadians();

    var a = 6377563.396, b = 6356256.909;      // Airy 1830 major & minor semi-axes
    var F0 = 0.9996012717;                     // NatGrid scale factor on central meridian
    var varphi0 = (49).toRadians(), lamda0 = (-2).toRadians();  // NatGrid true origin is 49°N,2°W
    var N0 = -100000, E0 = 400000;             // northing & easting of true origin, metres
    var e2 = 1 - (b*b)/(a*a);                  // eccentricity squared
    var n = (a-b)/(a+b), n2 = n*n, n3 = n*n*n; // n, n², n³

    var cosvarphi = Math.cos(varphi), sinvarphi = Math.sin(varphi);
    var upsilon = a*F0/Math.sqrt(1-e2*sinvarphi*sinvarphi);            // nu = transverse radius of curvature
    var rho = a*F0*(1-e2)/Math.pow(1-e2*sinvarphi*sinvarphi, 1.5); // rho = meridional radius of curvature
    var eta2 = upsilon/rho-1;                                    // eta = ?

    var Ma = (1 + n + (5/4)*n2 + (5/4)*n3) * (varphi-varphi0);
    var Mb = (3*n + 3*n*n + (21/8)*n3) * Math.sin(varphi-varphi0) * Math.cos(varphi+varphi0);
    var Mc = ((15/8)*n2 + (15/8)*n3) * Math.sin(2*(varphi-varphi0)) * Math.cos(2*(varphi+varphi0));
    var Md = (35/24)*n3 * Math.sin(3*(varphi-varphi0)) * Math.cos(3*(varphi+varphi0));
    var M = b * F0 * (Ma - Mb + Mc - Md);              // meridional arc

    var cos3varphi = cosvarphi*cosvarphi*cosvarphi;
    var cos5varphi = cos3varphi*cosvarphi*cosvarphi;
    var tan2varphi = Math.tan(varphi)*Math.tan(varphi);
    var tan4varphi = tan2varphi*tan2varphi;

    var I = M + N0;
    var II = (upsilon/2)*sinvarphi*cosvarphi;
    var III = (upsilon/24)*sinvarphi*cos3varphi*(5-tan2varphi+9*eta2);
    var IIIA = (upsilon/720)*sinvarphi*cos5varphi*(61-58*tan2varphi+tan4varphi);
    var IV = upsilon*cosvarphi;
    var V = (upsilon/6)*cos3varphi*(upsilon/rho-tan2varphi);
    var VI = (upsilon/120) * cos5varphi * (5 - 18*tan2varphi + tan4varphi + 14*eta2 - 58*tan2varphi*eta2);

    lamda = lamda-lamda0;
    var lamda2 = lamda*lamda, lamda3 = lamda2*lamda, lamda4 = lamda3*lamda, lamda5 = lamda4*lamda, lamda6 = lamda5*lamda;

    var N = I + II*lamda2 + III*lamda4 + IIIA*lamda6;
    var E = E0 + IV*lamda + V*lamda3 + VI*lamda5;

    return new OsGridRef(E, N);
};

/**
 * Ellipsoid parameters and datum parameters for transforming lat/lon coordinates between different
 * coordinate systems.
 *
 * @namespace
 */
var GeoParams = {};

/**
 * Ellipsoid parameters; major axis (a), minor axis (b), and flattening (f) for each ellipsoid.
 */
GeoParams.ellipsoid = {
    WGS84:        { a: 6378137,     b: 6356752.3142,   f: 1/298.257223563 },
    GRS80:        { a: 6378137,     b: 6356752.314140, f: 1/298.257222101 },
    Airy1830:     { a: 6377563.396, b: 6356256.909,    f: 1/299.3249646   },
    AiryModified: { a: 6377340.189, b: 6356034.448,    f: 1/299.32496     },
    Intl1924:     { a: 6378388.000, b: 6356911.946,    f: 1/297.0         },
    Bessel1841:   { a: 6377397.155, b: 6356078.963,    f: 1/299.152815351 }
};

/**
 * Datums; with associated *ellipsoid* and Helmert *transform* parameters to convert from WGS84
 * into given datum.
 */
GeoParams.datum = {
    WGS84: {
        ellipsoid: GeoParams.ellipsoid.WGS84,
        transform: { tx:    0.0,    ty:    0.0,     tz:    0.0,    // m
                     rx:    0.0,    ry:    0.0,     rz:    0.0,    // sec
                      s:    0.0 }                                  // ppm
    },
    OSGB36: { // www.ordnancesurvey.co.uk/docs/support/guide-coordinate-systems-great-britain.pdf
        ellipsoid: GeoParams.ellipsoid.Airy1830,
        transform: { tx: -446.448,  ty:  125.157,   tz: -542.060,  // m
                     rx:   -0.1502, ry:   -0.2470,  rz:   -0.8421, // sec
                      s:   20.4894 }                               // ppm
    },
    ED50: { // og.decc.gov.uk/en/olgs/cms/pons_and_cop/pons/pon4/pon4.aspx
        ellipsoid: GeoParams.ellipsoid.Intl1924,
        transform: { tx:   89.5,    ty:   93.8,     tz:  123.1,    // m
                     rx:    0.0,    ry:    0.0,     rz:    0.156,  // sec
                      s:   -1.2 }                                  // ppm
    },
    Irl1975: { // maps.osni.gov.uk/CMS_UserFiles/file/The_irish_grid.pdf
        ellipsoid: GeoParams.ellipsoid.AiryModified,
        transform: { tx: -482.530,  ty:  130.596,   tz: -564.557,  // m
                     rx:   -1.042,  ry:   -0.214,   rz:   -0.631,  // sec
                      s:   -8.150 }                                // ppm
    },
    TokyoJapan: { // www.geocachingtoolbox.com?page=datumEllipsoidDetails
        ellipsoid: GeoParams.ellipsoid.Bessel1841,
        transform: { tx:  148,      ty: -507,       tz: -685,      // m
                     rx:    0,      ry:    0,       rz:    0,      // sec
                      s:    0 }                                    // ppm
    }
};


/**
 * Creates lat/lon (polar) point with latitude & longitude values and height above ellipsoid, on a
 * specified datum.
 *
 * @classdesc Library of geodesy functions for operations on an ellipsoidal earth model.
 * @requires GeoParams
 * @requires Vector3d
 *
 * @constructor
 * @param {number}          lat - Geodetic latitude in degrees.
 * @param {number}          lon - Longitude in degrees.
 * @param {GeoParams.datum} [datum=WGS84] - Datum this point is defined within.
 * @param {number}          [height=0] - Height above ellipsoid, in metres.
 */
function LatLonE(lat, lon, datum, height) {
    if (typeof datum == 'undefined') datum = GeoParams.datum.WGS84;
    if (typeof height == 'undefined') height = 0;

    this.lat = Number(lat);
    this.lon = Number(lon);
    this.datum = datum;
    this.height = Number(height);
}


/**
 * Converts ‘this’ lat/lon coordinate to new coordinate system.
 *
 * @param   {GeoParams.datum} toDatum - Datum this coordinate is to be converted to.
 * @returns {LatLonE} This point converted to new datum.
 */
LatLonE.prototype.convertDatum = function(toDatum) {
    var oldLatLon = this;
    var transform;

    if (oldLatLon.datum == GeoParams.datum.WGS84) {
        // converting from WGS84
        transform = toDatum.transform;
    }
    if (toDatum == GeoParams.datum.WGS84) {
        // converting to WGS84; use inverse transform (don't overwrite original!)
        transform = {};
        for (var param in oldLatLon.datum.transform) {
            transform[param] = -oldLatLon.datum.transform[param];
        }
    }
    if (typeof transform == 'undefined') {
        // neither this.datum nor toDatum are WGS84: convert this to WGS84 first
        oldLatLon = this.convertDatum(GeoParams.datum.WGS84);
        transform = toDatum.transform;
    }

    // convert polar to cartesian
    var cartesian = oldLatLon.toCartesian();

    // apply transform
    cartesian = cartesian.applyTransform(transform);

    // convert cartesian to polar
    var newLatLon = cartesian.toLatLon(toDatum);

    return newLatLon;
};


/**
 * Converts ‘this’ point from polar (lat/lon) coordinates to cartesian (x/y/z) coordinates.
 *
 * @returns {Vector3d} Vector pointing to lat/lon point, with x, y, z in metres from earth centre.
 */
LatLonE.prototype.toCartesian = function() {
    var varphi = this.lat.toRadians(), lamda = this.lon.toRadians(), H = this.height;
    var a = this.datum.ellipsoid.a, b = this.datum.ellipsoid.b;

    var sinvarphi = Math.sin(varphi), cosvarphi = Math.cos(varphi);
    var sinlamda = Math.sin(lamda), coslamda = Math.cos(lamda);

    var eSq = (a*a - b*b) / (a*a);
    var upsilon = a / Math.sqrt(1 - eSq*sinvarphi*sinvarphi);

    var x = (upsilon+H) * cosvarphi * coslamda;
    var y = (upsilon+H) * cosvarphi * sinlamda;
    var z = ((1-eSq)*upsilon + H) * sinvarphi;

    var point = new Vector3d(x, y, z);

    return point;
};
/**
 * Creates a 3-d vector.
 *
 * The vector may be normalised, or use x/y/z values for eg height relative to the sphere or
 * ellipsoid, distance from earth centre, etc.
 *
 * @classdesc Tools for manipulating 3-d vectors, to support various latitude/longitude functions.
 *
 * @constructor
 * @param {number} x - X component of vector.
 * @param {number} y - Y component of vector.
 * @param {number} z - Z component of vector.
 */
function Vector3d(x, y, z) {
    this.x = Number(x);
    this.y = Number(y);
    this.z = Number(z);
}


/**
 * Returns the sum of ‘this’ vector and supplied vector.
 *
 * @param   {Vector3d} v - Vector to be added to this vector.
 * @returns {Vector3d} Vector representing sum of this and v.
 */
Vector3d.prototype.plus = function(v) {
    return new Vector3d(this.x+v.x, this.y+v.y, this.z+v.z);
};


/**
 * Returns the difference between ‘this’ vector and supplied vector.
 *
 * @param   {Vector3d} v - Vector to be subtracted from this vector.
 * @returns {Vector3d} This minus v.
 */
Vector3d.prototype.minus = function(v) {
    return new Vector3d(this.x-v.x, this.y-v.y, this.z-v.z);
};


/**
 * Returns the dot (scalar) product of ‘this’ vector and supplied vector.
 *
 * @param   {Vector3d} v - Vector to be dotted with this vector.
 * @returns {number} Dot product between ‘this’ and v.
 */
Vector3d.prototype.dot = function(v) {
    return this.x*v.x + this.y*v.y + this.z*v.z;
};


/**
 * Returns the cross (vector) product of ‘this’ vector and supplied vector.
 *
 * @param   {Vector3d} v - Vector to be crossed with this vector.
 * @returns {Vector3d} Cross product of ‘this’ and v.
 */
Vector3d.prototype.cross = function(v) {
    var x = this.y*v.z - this.z*v.y;
    var y = this.z*v.x - this.x*v.z;
    var z = this.x*v.y - this.y*v.x;

    return new Vector3d(x, y, z);
};


/**
 * Multiplies ‘this’ vector by a scalar value.
 *
 * @param   {number}   x - Factor to multiply this vector by.
 * @returns {Vector3d} Vector scaled by x.
 */
Vector3d.prototype.times = function(x) {
    return new Vector3d(this.x * x, this.y * x, this.z * x);
};


/**
 * Negates a vector to point in the opposite direction
 *
 * @returns {Vector3d} Negated vector.
 */
Vector3d.prototype.negate = function() {
    return new Vector3d(-this.x, -this.y, -this.z);
};


/**
 * Length (magnitude or norm) of ‘this’ vector
 *
 * @returns {number} Magnitude of this vector.
 */
Vector3d.prototype.length = function() {
    return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z);
};


/**
 * Normalizes a vector to its unit vector
 * – if the vector is already unit or is zero magnitude, this is a no-op.
 *
 * @returns {Vector3d} Normalised version of this vector.
 */
Vector3d.prototype.unit = function() {
    var norm = this.length();
    if (norm == 1) return this;
    if (norm === 0) return this;

    var x = this.x/norm;
    var y = this.y/norm;
    var z = this.z/norm;
    
    return new Vector3d(x, y, z);
};


/**
 * Calculates the angle between ‘this’ vector and supplied vector.
 *
 * @param   {Vector3d} v
 * @returns {number} Angle (in signed radians) between this vector and supplied vector.
 */
Vector3d.prototype.angleTo = function(v) {
    var sintheta = this.cross(v).length();
    var costheta = this.dot(v);

    return Math.atan2(sintheta, costheta);
};


/**
 * Rotates ‘this’ point around an axis by a specified angle.
 *
 * @param   {Vector3d} axis - The axis being rotated around.
 * @param   {number}   theta - The angle of rotation (in radians).
 * @returns {Vector3d} The rotated point.
 */
Vector3d.prototype.rotateAround = function(axis, theta) {
    // en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    // en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
    var p1 = this.unit();
    var p = [ p1.x, p1.y, p1.z ]; // the point being rotated
    var a = axis.unit();          // the axis being rotated around
    var s = Math.sin(theta);
    var c = Math.cos(theta);
    // quaternion-derived rotation matrix
    var q = [ [ a.x*a.x*(1-c) + c,     a.x*a.y*(1-c) - a.z*s, a.x*a.z*(1-c) + a.y*s],
        [ a.y*a.x*(1-c) + a.z*s, a.y*a.y*(1-c) + c,     a.y*a.z*(1-c) - a.x*s],
        [ a.z*a.x*(1-c) - a.y*s, a.z*a.y*(1-c) + a.x*s, a.z*a.z*(1-c) + c    ] ];
    // multiply q × p
    var qp = [0, 0, 0];
    for (var i=0; i<3; i++) {
        for (var j=0; j<3; j++) {
            qp[i] += q[i][j] * p[j];
        }
    }
    var p2 = new Vector3d(qp[0], qp[1], qp[2]);
    return p2;
    // qv en.wikipedia.org/wiki/Rodrigues'_rotation_formula...
};


/**
 * String representation of vector.
 *
 * @param   {number} [precision=3] - Number of decimal places to be used.
 * @returns {string} Vector represented as [x,y,z].
 */
Vector3d.prototype.toString = function(precision) {
    if (typeof precision == 'undefined') precision = 3;

    var p = Number(precision);

    var str = '[' + this.x.toFixed(p) + ',' + this.y.toFixed(p) + ',' + this.z.toFixed(p) + ']';

    return str;
};


/**
 * Converts ‘this’ point from cartesian (x/y/z) coordinates to polar (lat/lon) coordinates on
 * specified datum.
 *
 * @augments Vector3d
 * @param {GeoParams.datum.transform} datum - Datum to use when converting point.
 */
Vector3d.prototype.toLatLon = function(datum) {
    var x = this.x, y = this.y, z = this.z; 
    var varphi;

    var a = datum.ellipsoid.a, b = datum.ellipsoid.b;

    var eSq = (a*a - b*b) / (a*a);
    var p = Math.sqrt(x*x + y*y);
    varphi = Math.atan2(z, p*(1-eSq));

    var precision = 1 / a;  // 1m: Helmert transform cannot generally do better than a few metres
    do {
        var sinvarphi = Math.sin(varphi);
        var upsilon = a / Math.sqrt(1 - eSq*sinvarphi*sinvarphi);
        varphi = varphi;
        varphi = Math.atan2(z + eSq*upsilon*sinvarphi, p);
    } while (Math.abs(varphi-varphi) > precision);

    var lamda = Math.atan2(y, x);
    var H = p/Math.cos(varphi) - upsilon;

    var point = new LatLonE(varphi.toDegrees(), lamda.toDegrees(), datum, H);

    return point;
};

/**
 * Applies Helmert transform to ‘this’ point using transform parameters t.
 *
 * @private
 * @augments Vector3d
 * @param {GeoParams.datum.transform} t - Transform to apply to this point.
 */
Vector3d.prototype.applyTransform = function(t)   {
    var x1 = this.x, y1 = this.y, z1 = this.z;

    var tx = t.tx, ty = t.ty, tz = t.tz;
    var rx = (t.rx/3600).toRadians(); // normalise seconds to radians
    var ry = (t.ry/3600).toRadians(); // normalise seconds to radians
    var rz = (t.rz/3600).toRadians(); // normalise seconds to radians
    var s1 = t.s/1e6 + 1;             // normalise ppm to (s+1)

    // apply transform
    var x2 = tx + x1*s1 - y1*rz + z1*ry;
    var y2 = ty + x1*rz + y1*s1 - z1*rx;
    var z2 = tz - x1*ry + y1*rx + z1*s1;

    var point = new Vector3d(x2, y2, z2);

    return point;
};

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

/** Trims whitespace from string (q.v. blog.stevenlevithan.com/archives/faster-trim-javascript) */
if (typeof String.prototype.trim == 'undefined') {
    String.prototype.trim = function() {
        return this.replace(/^\s\s*/, '').replace(/\s\s*$/, '');
    };
}

/** Pads a number with sufficient leading zeros to make it w chars wide */
if (typeof String.prototype.padLz == 'undefined') {
    Number.prototype.padLz = function(w) {
        var n = this.toString();
        var l = n.length;
        for (var i=0; i<w-l; i++) n = '0' + n;
        return n;
    };
}

/** Extend Number object with method to convert numeric degrees to radians */
if (typeof Number.prototype.toRadians == 'undefined') {
    Number.prototype.toRadians = function() { 
      return this * Math.PI / 180; 
    };
}

/** Extend Number object with method to convert radians to numeric (signed) degrees */
if (typeof Number.prototype.toDegrees == 'undefined') {
    Number.prototype.toDegrees = function() { 
      return this * 180 / Math.PI; 
    };
}


/**
* Espruino GPS
* 
*/
function onInit() {
	SPI1.setup({ sck:B3, mosi:B5 });
	var g = require("PCD8544").connect(SPI1,B6,B7,B8, function() {
		g.clear();
		g.drawString("GPS!",0,0);
		g.drawLine(0,5,84,5);
		g.flip();
		console.log('screen initialised');
	});

    Serial4.setup(9600,{tx:C10,rx:C11});
	var gps = require("GPS").connect(Serial4, function(data) {
		//console.log(data);
		g.clear();
		g.setFontBitmap();
		g.drawString("Sat:" + data.satellites + "(" + nelmt.length + ")", 0,0);

		if (data.fix == 1) {
			var os = OsGridRef.latLongToOsGrid(data);
			//console.log(data);
			nelmt.push(os.northing);
			//console.log(nelmt.length);
			if (nelmt.length == 10) {
				navg = '';
				nsum = 0;
				for(var i = 0; i < nelmt.length; i++) {
					nsum += parseInt(nelmt[i], 10);
				}
				navg = Math.round(nsum/nelmt.length);
				while(nelmt.length > 0) {
					nelmt.pop();
				}
			}
			eelmt.push(os.easting);
			if (eelmt.length == 10) {
				eavg = '';
				esum = 0;
				for(var t = 0; t < eelmt.length; t++) {
					esum += parseInt(eelmt[t], 10);
				}
				eavg = Math.round(esum/eelmt.length);
				while(eelmt.length > 0) {
					eelmt.pop();
				}
			}
			g.drawString("Alt:" + data.altitude, 35,0);
			g.drawLine(0,5,84,5);
			g.setFontVector(12);
			g.drawString("N " + navg, 0,10);
			g.drawString("E " + eavg, 0,29);
			//console.log('northing ' + os.northing + ' easting ' +os.easting);  
		}else if (data.fix === 0){
			g.drawString("Sat Fix...", 35,0);
            g.drawLine(0,5,84,5);
			if (navg > 0) {
				g.setFontVector(12);
				g.drawString("N " + navg, 0,10);
				g.drawString("E " + eavg, 0,29);
			}

		}
		g.flip();
	});	
}
onInit();
