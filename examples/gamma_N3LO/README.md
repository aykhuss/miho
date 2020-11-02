# Off-shell photon productin at N3LO

## Prediction using the ABC model

### 1. vanilla example
run
```
miho abc
```
and enter the cross section numbers
```
7.110012e+00 8.192043e+00 8.154210e+00 7.956180e+00 
```
and the code will spit out three numbers (median, lower & upper bound of the 68% DoB interval):
```
7.941417556640626 7.88936051953125 7.995028535156253
```
To embed it into scripts, you might want to pipe:
```
echo "7.110012e+00 8.192043e+00 8.154210e+00 7.956180e+00" | miho abc
```


### 2. using file as input
If you have a file that contains the cross section numbers (single space-separated line), you can also run
```
miho -f data/N3LO_scl1.0-1.0.dat abc
```


### 3. setting model parameters
The ABC model has 3 parameters that fix its behaviour. 
You can set those via options *after* the model subcommand (help display: `miho <subcommand> -h`)
```
miho abc --omega 2 --xi 0.5 --epsilon 0.3
```
or using an input file
```
miho -f data/N3LO_scl1.0-1.0.dat abc --omega 2 --xi 0.5 --epsilon 0.3
```


### 4. controling the output
By default we get the median and the 68% DoB interval but using the format option we can control what we want (full list of template variables in main README.md) and how it's displayed:
```
miho --format "{mean} +/- {stdev} -- {median} [{dob95_low}, {dob95_upp}]" -f data/N3LO_scl1.0-1.0.dat abc
```
should give
```
7.9315555750799644 +/- 0.25853020905374113 -- 7.941417556640626 [7.682686312500001, 8.2296736875]
```


### 5. get the PDF
If you want to have the probability distribution (list of x-y value pairs), add the `--pdf` flag and the output of the nodes of the importance sampling will be written out.
If you want more nodes to be outputted, set the `--nmax <int>` and `--accuracy <float>` and request higher precision. 
```
miho --pdf -f data/N3LO_scl1.0-1.0.dat abc
```
If you request the PDF, the format will be re-set. 
In case you want more information to be displayed, set the `--format` option as described in pt. 4.


### 6. one-dimensional scale marginalisation
Let's try to marginalise over 1 dimension using the "Alekas prescription" and the Gauss-Legendre method.
We do this by passing the `F` factors to the `--scl1d` option and setting the `--gl` flag.
You will now be asked to provide numbers for the different scales separately:
```
miho --scl1d 0.5 1.0 2.0 --gl abc
# ABCNumericModel
entered app_scl1d
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
Clearly, this becomes quite cumbersome so you can provide an input file like this:
```
cat data/N3LO_3scl.dat
0.5  6.187361e+00 7.956126e+00 8.106391e+00 7.893847e+00
1.0  7.110012e+00 8.192043e+00 8.154210e+00 7.956180e+00
2.0  7.966379e+00 8.409684e+00 8.231636e+00 8.020915e+00
```
that contains the scale factor `F` at the beginning of the line followed by the cross-section numbers.
Then running
```
miho -f data/N3LO_3scl.dat --scl1d 0. --gl abc
# ABCNumericModel
7.940640585937502 7.816325273437501 8.061848015625003
```
(note: I need to provide some number(s) to `--scl1d` for it to parse correctly. Will be made more elegant )


### 7. two-dimensional scale marginalisation
Now we need to provide scale pairs `--scl2d [1.,1.] [0.5,1.] ...`.
In case you use input files, the format is:
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
and the output is
```
miho -f data/N3LO_9scl.dat --scl2d [0.,0.] --gl abc
# ABCNumericModel
7.921993289062502 7.732412437500001 8.117789906250001
```




