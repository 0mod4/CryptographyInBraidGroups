Masterthesis in mathematics, cryptography
=========================================
Cryptography in braid groups
----------------------------
0mod4

Group based cryptography describes a number of cryptographic methods based on group theoretic problems. One of these methods is the Anshel-Anshel-Goldfeld protocol (AAG) which is based on the generalized conjugacy search problem in braid groups. This thesis considers length based attacks on the AAG protocol as treated in 'Length-based conjugacy search in the braid group' by D. Garber, S. Kaplan, M. Teicher, B. Tsaban and U. Vishne.
 
To do so, first a short overview on group based Cryptography, braid groups and the AAG-protocol is given. The book 'Group-based Cryptography' by A. Myasnikov, V. Shpilrain and A. Ushakov was inspiration and help for these parts.

Since the algorithm used to obtain random elements in a braid group used in the paper does not generate uniformly distributed random elements within a given range of lengths, some other methods to generate random elements are considered. 
This is done to see how the test might perform in the 'real world' and to see if uniformly distributed random input makes it more or less likely/fast to break the AAG protocol.
In particular the algorithm to generate uniform elements of given lengths as described by V. Gebhardt and J. Gonzáles-Meneses in 'Generating random braids' is implemented and used throughout all experiments.

The executed experiments are inspired by 'Length-based conjugacy search in the braid group', additional tests are performed, for example a closer look at the runtime of length based attacks. 
Gained results are then discussed in the last section. The optimal attack under the tested ones is found, the security of the AAG-scheme under this attack is considered and the results from the paper are compared to the ones
using the random generator from 'Generating random braids'.

All code needed for this thesis can be found in this repository
