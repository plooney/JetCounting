### JetCounting: Image Analysis Technique for Measurement of 3D Spiral Artery Jets using 3D power Doppler Ultrasound

__JetCounting__ is a tool for the analysis of 3D power Doppler ultrasound image data. Expected inputs are as follows. 

```bash
JetCounting inputBModeVolume inputPowerDopplerVolume inputSegVolume inputDTVolume dist outputDir [axis] [interface]
```
Runs the JetCounting main Application. Using 4 images that are pre-processed the application will extract the poly data meshes at defined distance in millimeters from the utero-placenta interface and then output the 3D image files and meshes to an output directory.

Choice of 3D axis to extract and orientation can be user defined by the optional arguments.

**this program assumes all volumes are in the same resolution and spacing, in MetaImage format. NB Resampling may be required. **

- inputBModeVolume path to the 3D bmode volume
- inputPowerDopplerVolume path to the 3D PD volume
- inputSegVolume path to the 3D binary segmentation volume
- inputDTVolume path to the 3D Signed Euclidean Distance Transform volume
- dist floating point argument of distance from the UPI to measure
- outputDir output directory to store volumes and meshes for analysis
- axis optional integer parameter of axis (0,1,2) in which extraction should occur
- interface optional integer parameter of negative or positive normal (0,1) extraction of the UPI from the total mesh so occur.

## Requirements/Dependencies
This program requires VTK 7.0 and ITK 4.8.2. It has been successfully built on Windows 7 using Visual Studio 12.