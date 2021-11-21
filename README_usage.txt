
This app will take a280 absorbance values returned from the FDB UV-Vis plate reader protocol and return a heatmap and a CSV in mg/ml

It will use protein metadata from the proteinsequences.csv to calculate:
- MW in KD
- Extinction Coefficent reduced

It works on the basis that 20 ul of elution is mixed with 80ul of HT buffer and results in a well volume of 100ul with a path length of 0.2857 cm.

1. Put the raw data file in the data directory
2. YOU MUST ENTER THE PROTEIN NAMES AS PER THE PROTEINSEQUENCES CSV IN THE ORDER THAT THEY APPEAR IN THE WELLS A THROUGH H.
