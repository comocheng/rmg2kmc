# Code to convert RMG mechanism files into KMC input files

## Planned KMC Codes to Support
- Monte Coffee
- Zacros
- SPPARKS

## TODO
- Cantera mechanism reader
- Chemkin mechanism reader
- Monte Coffee writer
- Zacros writer
- identify favored binding sites - use RMG species dictionary to identify which atom is bound to the surface?
- specify multiple surface facets (to start, assume fcc(111) facet)

## Usage
1. Create a directory with your mechanism file
2. `python rmg2kmc.py input_file.yaml`

