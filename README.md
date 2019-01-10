# nf2ff_transformation

This is a set of Matlab scripts and functions to transform near-field 
antenna meusurements into far-field.

## Transformation Scripts
For each scanning setup (planar,cylindrical,spherical) there is a main script that handles the transformation of the near-field measurements.
To run one of them, adjust the path to the measurement data files inside the first section of the script.
Then select whether you want to use the FFT implementation or the so-called matrix method (by-hand-calculation). 
The former being the faster of course.
Afterwards, several different plots can be created. 
Either pattern cuts at different angles, the accumulated pattern error, the error in the half-power beamwidth or the error in the first side-lobe level.
In case of the script for the [probe position noise analysis](nf2ff_planar_noisy.m) the above mentioned errors are plotted over the standard devitiation of the positioning noise.


## Additional Functions
For visiualization several additional functions have been implemented.
Using [plotPlanarNFData](plotPlanarNFData.m) and [plotCylindricalNFData](plotCylindricalNFData.m) the data that is read from the measurement files can be plotted.

Using the [TransformationResults](TransformationResults.m) function, the results of a specific scan can be visualized in 3D either in cartesian or spherical coordinates. 
Additionally the difference between the two pattern is plotted.
