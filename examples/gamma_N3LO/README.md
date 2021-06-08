# Off-shell photon productin at N3LO

## Prediction using the abc model

### 1. vanilla example
Execute the program
```
miho abc
```
and enter the cross section numbers in the prompt
```
7.110012e+00 8.192043e+00 8.154210e+00 7.956180e+00 
```
The code will display three numbers that correspond to the median and the lower & upper bound of the 68% DoB interval:
```
7.941417556640626 7.88936051953125 7.995028535156253
```
Embedding the tool into existing workflows can be accomplished by a simple pipe:
```
echo "7.110012e+00 8.192043e+00 8.154210e+00 7.956180e+00" | miho abc
```


### 2. using a file as input
The cross section numbers can also be read in from a file that contains the cross section numbers (single space-separated line):
```
miho -f data/N3LO_scl1.0-1.0.dat abc
```


### 3. setting the model parameters
The abc model has three parameters that fix its properties. 
Those can be set via options *after* the model subcommand (help display: `miho <subcommand> -h`):
```
miho abc --omega 2 --xi 0.5 --epsilon 0.3
```
or in the case of using an input file (before the model specification):
```
miho -f data/N3LO_scl1.0-1.0.dat abc --omega 2 --xi 0.5 --epsilon 0.3
```


### 4. controlling the output
By default, the median and the 68% DoB intervals are returned but using the format specification this can be easily adjusted (full list of template variables in the main README.md):
```
miho --format "{mean} +/- {stdev} -- {median} [{dob95_low}, {dob95_upp}]" -f data/N3LO_scl1.0-1.0.dat abc
```
should give
```
7.9315555750799644 +/- 0.25853020905374113 -- 7.941417556640626 [7.682686312500001, 8.2296736875]
```


### 5. the probability distribution
The probability distribution can be obtained by adding the flag `--pdf`. This will output the nodes of the importance sampling as a list of x-y value pairs.
In case more nodes are desired, the flags `--nmax <int>` and `--accuracy <float>` allow to adjust the precision targets: 
```
miho --pdf -f data/N3LO_scl1.0-1.0.dat abc
```
Note that in the case the PDF is requested, the format will be re-set. 
In case additional information is needed, the `--format` option as described in pt. 4 above must be used.


### 6. one-dimensional scale marginalisation
We consider the case with a one-dimensional scale dependence that we treat in the "scale averaging" procedure using Gauss-Legendre quadrature.
This is done by passing the `F` factors to the `--scl1d` option and setting the `--gl` flag.
The prompt will ask to provide numbers for the different scales separately:
```
miho --scl1d 0.5 1.0 2.0 --gl abc
# ABCNumericModel
# Enter XS[0.5×μ₀] values @ LO NLO ... separated by spaces:
6.187361e+00 7.956126e+00 8.106391e+00 7.893847e+00
# add_model: _n_orders = 0
# Enter XS[1.0×μ₀] values @ LO NLO ... separated by spaces:
7.110012e+00 8.192043e+00 8.154210e+00 7.956180e+00
# add_model: _n_orders = 4
# Enter XS[2.0×μ₀] values @ LO NLO ... separated by spaces:
7.966379e+00 8.409684e+00 8.231636e+00 8.020915e+00
# add_model: _n_orders = 4
7.940640585937502 7.816325273437501 8.061848015625003
```
The equivalent input file
```
cat data/N3LO_3scl.dat
0.5  6.187361e+00 7.956126e+00 8.106391e+00 7.893847e+00
1.0  7.110012e+00 8.192043e+00 8.154210e+00 7.956180e+00
2.0  7.966379e+00 8.409684e+00 8.231636e+00 8.020915e+00
```
contains the scale factor `F` at the beginning of the line followed by the cross-section numbers.
```
miho -f data/N3LO_3scl.dat --scl1d 0. --gl abc
# ABCNumericModel
7.940640585937502 7.816325273437501 8.061848015625003
```
(note: some number(s) must be given to `--scl1d` for it to parse correctly but that are ignored)


### 7. two-dimensional scale marginalisation
The two-dimensional case proceeds by specifying pairs of scale multipliers `--scl2d [1.,1.] [0.5,1.] ...`.
The input file format is a generalisation of the one-dimensional case:
```
cat data/N3LO_9scl.dat
0.5 0.5  6.187361e+00 7.956126e+00 8.106391e+00 7.893847e+00
0.5 1.0  6.187361e+00 7.766542e+00 8.046989e+00 7.944202e+00
0.5 2.0  6.187361e+00 7.614409e+00 7.974627e+00 7.958674e+00
1.0 0.5  7.110012e+00 8.321943e+00 8.136375e+00 7.877205e+00
1.0 1.0  7.110012e+00 8.192043e+00 8.154210e+00 7.956180e+00
1.0 2.0  7.110012e+00 8.087803e+00 8.146807e+00 8.007180e+00
2.0 0.5  7.966379e+00 8.515837e+00 8.117082e+00 7.862292e+00
2.0 1.0  7.966379e+00 8.456943e+00 8.188999e+00 7.950480e+00
2.0 2.0  7.966379e+00 8.409684e+00 8.231636e+00 8.020915e+00
```
that should generate the output similar to
```
miho -f data/N3LO_9scl.dat --scl2d [0.,0.] --gl abc
# ABCNumericModel
7.921993289062502 7.732412437500001 8.117789906250001
```




