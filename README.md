# IIS Cutting generation algorithm for binary chance-constrained problems

This is the code used in the paper
 - Canessa, G., Gallego, J.A., Ntaimo, L. et al. An algorithm for binary linear chance-constrained problems using IIS. Comput Optim Appl 72, 589–608 (2019). https://doi.org/10.1007/s10589-018-00055-9

It will not be mantained, I just uploaded it for legacy purposes. However, feel free to ask.

Usage:
- Setting up the necessary packages
```shell
pip install -r requirements.txt
```

- Running the code
```shell
python bab.py model beta cuts
```
*model* is the name of the SMPS format file to be used (don't add the .cor/sto/tim extension).
*beta* is the value for the chance-constrained probability.
*cuts* if 1 then we will introduce the cutting scheme from the paper, 0 it will be a normal Gurobi run.