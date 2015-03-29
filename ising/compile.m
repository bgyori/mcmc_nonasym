fprintf('Compiling Curie-Weiss\n');
mex -output isingCW -O -IisingCW -Imc -I. isingCW/magIsingCW.cpp
fprintf('Compiling Curie-Weiss Metropolis\n');
mex -output isingCWMetropolis -O -DMETROPOLIS -IisingCW -Imc -I. isingCW/magIsingCW.cpp
fprintf('Compiling Curie-Weiss sign of magnetization\n');
mex -output isingCWMagsign -O -DMAGSIGN -IisingCW -Imc -I. isingCW/magIsingCW.cpp
fprintf('Compiling Ising 1D\n');
mex -output ising1D -O -Iising1D -Imc -I. ising1D/magIsing1D.cpp
fprintf('Compiling Ising 1D systematic scan\n');
mex -output ising1Dsyst -O -DSWEEP -Iising1D -Imc -I. ising1D/magIsing1D.cpp
fprintf('Compiling Ising 2D\n');
mex -output ising2D -O -Iising2D -Imc -I. ising2D/magIsing2D.cpp

