# collimator
This is a collection of scripts to design 3D printed radial ("Soller slit") collimators.

## Prerequisites
- Python (all components should now work with Python 3)
  - numpy
  - scipy
  - vapory (to allow raytracing with POV-Ray in Python)
- [FreeCAD](https://www.freecad.org/) 0.18
- [POV-Ray](https://www.poyray.org/) 3.7.0.8
- [stltools](https://github.com/rsmith-nl/stltools)

## Steps to design a Soller collimator
### 1. Design your grid
Edit the script `designer/soller.py` to meet your requirements. Most parameters are as in the earier Mathematica vesion of the script.
  - `format` should usually be `json`
    - `scad` is retained as an option to produce output compatible with [OpenSCAD](https://openscad.org/), which can allow a quick preview of the model without having to render
  - `part` should be `0` unless you want to generate only a subset of the grid walls (this can be used as part of the process to facilitate blob-free printing)

Run the script to produce the `RadialCollimatorBox` and `RadialCollimatorGrid` files.
  
### 2. Create the 3D model
  1. Open `Collimator.FCMacro` in FreeCAD.
      - Ensure the paths to your .json files are correct
      - Edit the relevant lines to set the desired wall thicknesses
  2. Create a new empty document
  3. Execute the macro in the editor
    
### 3. Test the model using raytracing
  1. Use `stl2pov` from stltools to produce a `.inc` file
  2. Use the `raytracing/povray.py` script, in particular its `scene` function to render what an image would look like on an area detector; examples can be found in the script. It is recommended to use an interactive Python shell while working with this
     - Always be sure to include the `auto_camera_angle=False` option, otherwise the view will not be as you intended

### 4. Print the grid
Your mileage may vary, depending on the material used and the print options selected.

## Further reading
- [A novel 3D printed radial collimator for x-ray diffraction](https://aip.scitation.org/doi/suppl/10.1063/1.5063520) by S. Kowarik, L. Bogula, S. Boitano, F. Carlà, H. Pithan, P. Schäfer, H. Wilming, A. Zykov and L. Pithan 

