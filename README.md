Protein Projected Area Analysis Tool
Overview
This Python script calculates the projected surface area of a protein onto a 2D plane after dual-axis rotation in 3D space, generating a matrix data file of rotation angles versus projected areas. The tool is designed for structural biology and molecular modeling applications, useful for analyzing how protein exposure changes from different viewing angles.
Features
• Extracts all alpha-carbon (CA) atom coordinates from a PDB file
•	Rotates the protein model around X and Y axes (0° to 180° with 0.5° increments)
•	Projects rotated 3D coordinates onto the XY plane
•	Calculates 2D projected area using a grid-based method (considering atom radius)
•	Outputs a 360×360 matrix file containing projected areas for all angle combinations
Requirements
•	Python 3.6 
•	NumPy
Installation
bash
pip install numpy
Usage
Basic Usage
bash
python protein_projection.py protein.pdb
Command Line Options
text
positional arguments:
  input_pdb             Input PDB file path

optional arguments:
  -o, --output OUTPUT   Output filename (default: area_matrix.dat)
  --step STEP           Grid step size in Ångströms (default: 1.0)
  --radius RADIUS       Atom radius in Ångströms (default: 1.7)
  -h, --help            Show help message
Examples
1.	Default settings:
bash
python protein_projection.py 1xyz.pdb
2.	Custom output file and grid precision:
bash
python protein_projection.py 1xyz.pdb -o area_results.txt --step 0.5
3.	Custom atom radius:
bash
python protein_projection.py protein.pdb --radius 2.0 --step 1.0
Output Format
The script generates a tab-separated matrix file where:
•	Rows correspond to X-axis rotation angles (0° to 179.5°, 0.5° increments)
•	Columns correspond to Y-axis rotation angles (0° to 179.5°, 0.5° increments)
•	Each cell contains the projected area in Å² for that specific orientation
Example header format:
text
Y000.0 Y000.5 Y001.0 ... Y179.5
[area values...]
Algorithm Details
1.	Coordinate Extraction: Reads CA atom coordinates from PDB ATOM records
2.	Rotation: Applies 3D rotation matrices around X and Y axes
3.	Projection: Drops Z-coordinate to project onto XY plane
4.	Area Calculation: Uses a grid-based method where each atom covers a circular area with specified radius
o	Creates a grid covering all projected atoms
o	Marks grid points within atom radius as covered
o	Sums covered area as: covered_grid_points × (step_size²)
Performance Considerations
•	The script processes 360 × 360 = 129,600 angle combinations
•	Each combination requires rotation, projection, and area calculation
•	Computation time scales with:
o	Number of CA atoms in the protein
o	Grid resolution (smaller step = more computation)
o	Protein size (larger proteins = larger grids)
Troubleshooting
1.	"No CA atoms found" error:
o	Ensure PDB file contains ATOM records with CA atoms
o	Check file format and integrity
2.	Large memory usage:
o	Reduce grid resolution with larger --step values
o	Consider processing smaller proteins
3.	Long computation time:
o	Increase grid step size
o	Process only a subset of angles by modifying the script
Citation
If you use this tool in your research, please cite appropriately based on your application context.
License
This tool is provided for academic and research use. Please contact the author for commercial applications.
Contact
For questions, suggestions, or bug reports, please refer to the original source repository or contact the developer directly.

