This folder contains the source code of the implementations experimented in the submission 16 of the Arith 27 conference:

Each subfolder corresponds to one of the experimentation described in the submission.

- TABLE1: comparison between the gf2x_mul and our KaratRecMult;
- TABLE2:
  - Toom3 : comparison between the gf2x_mul and our Toom3Mult;
  - Toom3Rec : comparison between the gf2x_mul and our Toom3RecMult (1-recursive version);
- TABLE3: comparison between the gf2x_mul and our Toom3Mult and Toom3RecMult for 47647 and 70853 HQC sizes;
- TABLE4: comparison between the HQC implemented sparse-dense multiplication and our FastConvolution (with and without randomization);
- TABLE5: comparison between the HQC patched with the gf2x_mul and with our Toom3Mult and Toom3RecMult;
- TABLE6: comparison with the NIST submitted optimzed implementation and our patched versions with FastConvolution multiplication (with and without randomization)

*Warnings:*
These source codes are designed for compilation and execution on x86-64 platforms, with AVX, AVX2 and PCLMUL extensions.

To run the experimentations of TABLE1, TABLE2, TABLE3 and TABLE5, you may need to install the gf2x library

-> https://gforge.inria.fr/frs/?group_id=1874

The experimentations of the submission have been made with gf2x 1.2 version of the 03/07/2017.
