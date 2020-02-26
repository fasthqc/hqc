- sparse_dense_mul is the multiplication of the original HQC source code;
- fastConvolutionMult is the multiplication of the Arith2020 submission 16.

This folder contains two subfolders:

- FC: source code for fastConvolutionMult vs sparse_dense_mul;
- RandFC: the same except the fastConvolutionMult which is the randomized version.

The procedure is the same for both subfolders.
To launch individual tests, please type:

    $ make -B SIZE=XX OMEGA=YY
    
The values for SIZE must be in the set {24677,43669,46747,63587,67699,70853}, otherwise, you may have bad results.

To launch all the tests, please type:

    $ make test
    
The results of the performances might vary between global tests and individual tests, and the best results are obtained with the individiual test procedure.
