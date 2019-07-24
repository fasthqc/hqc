AVX2 Implementation of HQC (Last update 24/07/2019).

We propose in this repository a modification of the optimized HQC source code (HQC NIST round 2 submission, 10/04/2019, see Aguilar-Melchioret al.[1]) in order to improve the performances.  This modification deals with the vectorization of two functions:  the *sparse-dense polynomial*  multiplication and  the  *syndrome generation*. More details can be found in the file *report.pdf*. 

The folders *Optimized_Implementation* and *Reference_Implementation 2* constain the source code of the original NIST proposal.
**Our imrovement is in the *Fast_implementation* folder.**

**The speed-up of our patched HQC implementation is about 34 to 46 % in comparison to the current optimized NIST proposal**

*[1]  Carlos Aguilar-Melchior, Nicolas Aragon, Slim Bettaieb, Loïc Bidoux, Olivier Blazy, Jean-Christophe Deneuville,  Philippe Gaborit,  Edoardo Persichetti,  and Gilles Zémor.   Hamming Quasi-Cyclic (HQC).  InNIST Post-Quantum Cryptography submissions, round 2. NIST,2019.*
