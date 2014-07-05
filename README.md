espruino-gps
============

<h2>Convert between GPS WGS84 to OSGB36 OS grid reference</h2>

<p>Ublox NEO6MV2 GPS + PCD8544 LCD DRIVER (NOKIA 5110)</p>

The GPS example at [http://www.espruino.com/Pocket+Walking+GPS] for UK OS grid co-ordinates doesn't take into consideration 
ellipsoid parameters and the differences in transforming lat/lon coordinates between different coordinate systems. The Ublox 
NEO6MV2 GPS will take down GPS WGS84 co-ordinates and must be converted to OSGB36 with an additional Helmert transform 
otherwise the OS grid reference will be out by upto 120m.

Thanks to Chris Veness for his Geodesy tools of an ellipsoidal earth model at [http://www.movable-type.co.uk/]

<h3>The outcome</h3>
1. Take 10 GPS readings
2. Average
3. Convert between GPS WGS84 to OSGB36
4. Display OS grid reference
