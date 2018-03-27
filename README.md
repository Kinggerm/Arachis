# Arachis

## Introduction

Arachis is a Python library for analyzing genome rearrangements. It allow users to reconstruct ancestral genome 
gene orders and infer pairwise genome differences or events.


## Algorithms & Features

The algorithm for reconstructing ancestral genome gene orders implemented in the script file `run_pypmag.py`
is derived from the ancestral gene order reconstruction module of <a href="#PMAG">`PMAG+`</a>, with modifications:

1. Circular and gap-containing genomes is allowed as inputs. See modification on GRIMM below.
2. Equipped with python multiprocessing.
3. More flexible in input both tree and data format.

This library defines a new version of [the classic GRIMM format](http://grimm.ucsd.edu/GRIMM/grimm_instr.html) 
with following modifications:
1. Blocks could be named with letters in "`-.0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ_abcdefghijklmnopqrstuvwxyz|~[]`". 
"`-`" still means a reverse direction when it appears at the first letter.
2. The "`*`" block in a sequence means a gap. `1 2 3 * 4 5 * 6 $` is equivalent to `1 2 3 * -5 -4 * 6 $`.
3. A sequence line without the "`$`" block in the end stands for a circular chromosome. A circular sequence
`a b c d e` is equivalent to `b c d e a`, also equivalent to `-e -d -c -b -a`. But in this case only 
one chromosome per sample is allowed, and multiple lines without "$" would be regarded as one single chromosome 
written interleaved. This design is due to the limitation of applying to tsp solver.

The functions in Arachis for inferring pairwise genome differences or events are still at <b>infant</b> stage. 
If you find any bugs or something to improve, please contact 
[jinjianjun@mail.kib.ac.cn](mailto:jinjianjun@mail.kib.ac.cn). New contributors are welcome! 
Also, users have to bear in mind that do not test data with too many breakpoints (like 10+).
Currently, the function SignedPermutation.inversion_event_from utilizes an <b>exhausted</b> 
scheme searching for one best solution. Currently, I'm using it to play with highly rearranged plastome data of legumes.
It's worth trying more small permutations, like some plant mitochondrial data.

## Installation

Download Arachis and install Arachis with:

    $ git clone "https://github.com/Kinggerm/ARACHIS"
    $ cd ARACHIS
    $ python setup.py install

To further use `run_pypmag.py` to reconstruct ancestral genome gene order, you have to install following dependencies:

* <b>DendroPy</b> The tree parser in Arachis. Get it [here](https://www.dendropy.org/#installing).
* <b>RAxML</b> The reconstruction engine in the algorithm of PMAG+. 
The single thread version is preferred. Get it [here](https://github.com/stamatak/standard-RAxML).
* <b>Concorde</b> The TSP (Traveling Salesman Problem) Solver in the algorithm of PMAG+. 
Get it [here](http://www.math.uwaterloo.ca/tsp/concorde/downloads/downloads.htm).

## Example

* To check whether two circular permutations, `-e -d -c -b -a` and `a b c d e`, are equivalent:

    ```py
        # open python shell
        >>> from arachis.genomeClass import Chromosome
        >>> seq1 = Chromosome("-e -d -c -b -a")
        >>> seq2 = Chromosome("a b c d e")
        >>> seq1 == seq2
    ```
     <pre>    <font color=grey>True</font></pre>

* If you want to see how many flip-flop configurations (isomers) could be induced by several groups of inverted repeats, 
or in another similar case, to see how many reasonable paths are there in a complicated assembly graph with repeats that
could not be unfolded by short seq-library, try this:
 
     ```py
        # open python shell
        >>> from arachis.genomeClass import Chromosome
        >>> Picea = Chromosome("1 2 12 14 13 2 3 4 10 8 15 14 11 4 5 6 7 8 9 -6")
        >>> isomers, changes = Picea.get_isomers()
        >>> print(len(isomers))
     ```
     ```
        15
     ```
     ```py
        >>> for isomer in isomers:
                print(isomer)
     ```
     ```
        1 2 12 14 13 2 3 4 10 8 15 14 11 4 5 6 -9 -8 -7 -6
        1 2 12 14 13 2 3 4 10 8 9 -6 -5 -4 -11 -14 -15 -8 -7 -6
        1 2 12 14 11 4 5 6 -9 -8 -10 -4 -3 -2 -13 -14 -15 -8 -7 -6
        1 2 12 14 13 2 3 4 5 6 -9 -8 -10 -4 -11 -14 -15 -8 -7 -6
        1 2 3 4 10 8 9 -6 -5 -4 -11 -14 -12 -2 -13 -14 -15 -8 -7 -6
        1 2 12 14 11 4 10 8 9 -6 -5 -4 -3 -2 -13 -14 -15 -8 -7 -6
        1 2 12 14 11 4 5 6 7 8 15 14 13 2 3 4 10 8 9 -6
        1 2 12 14 13 2 3 4 5 6 7 8 15 14 11 4 10 8 9 -6
        1 2 3 4 5 6 -9 -8 -10 -4 -11 -14 -12 -2 -13 -14 -15 -8 -7 -6
        1 2 3 4 10 8 15 14 13 2 12 14 11 4 5 6 -9 -8 -7 -6
        1 2 12 14 11 4 10 8 15 14 13 2 3 4 5 6 -9 -8 -7 -6
        1 2 3 4 5 6 7 8 15 14 13 2 12 14 11 4 10 8 9 -6
        1 2 3 4 10 8 15 14 13 2 12 14 11 4 5 6 7 8 9 -6
        1 2 12 14 11 4 10 8 15 14 13 2 3 4 5 6 7 8 9 -6
    ```

* Run `run_pypmag.py` to reconstruct ancestral genome gene order of test data:

        run_pypmag.py -d test/test_1_grimm.txt -t test/test_1_rooted.tre -o test/test_1_output --seed 12345
    
* To see parsimonious events along the branch from `A1` to `sp2` in above test_1 running results:

     ```py
        # open python shell
        >>> from arachis.genomeClass import GenomeList
        >>> extant_samples = GenomeList("test/test_1_grimm.txt")
        >>> sp2 = extant_samples["sp2"].chromosomes()[0]
        >>> ancestors = GenomeList("test/test_1_output/OutputGeneOrder")
        >>> A1 = ancestors["A1"].chromosomes()[0]
        >>> events = sp2.event_from(A1)
     ```
     ```
                Breakpoints: 2
                    Round 1: inherited combinations: 1; inversion sites:  2; time: 0.0002s; memory: 0.01G
                Inversions: 1 + 0(iso)
                Total inversion time: 0.0006s
     ```

## Citation

If you use Arachis in your research, you could cite Arachis as:
* Jian-Jun Jin. 2018. Arachis: a Python library for analysing genome rearrangements. https://www.github.org/Kinggerm/ARACHIS

If you use `run_pypmag.py`, please cite following papers:
<div id="PMAG"></div>

* <b>PMAG+</b> Hu, F., Zhou, J., Zhou, L., & Tang, J. 2014. Probabilistic Reconstruction of Ancestral Gene Orders with Insertions and Deletions. IEEE/ACM Transactions on Computational Biology and Bioinformatics, 11(4), 667–672. <http://doi.org/10.1109/TCBB.2014.2309602>
* <b>RAxML</b> Stamatakis, A. 2014. RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics, 30(9), 1312–1313. <http://doi.org/10.1093/bioinformatics/btu033>
* <b>Concorde</b> D. Applegate, R. Bixby, V. Chvatal, & W. Cook. 2011. “Concorde TSP Solver,” http://www.math.uwaterloo.ca/tsp/concorde/
* <b>DendroPy</b> Sukumaran, J., & Holder, M. T. 2010. DendroPy: a Python library for phylogenetic computing. Bioinformatics, 26(12), 1569–1571. <http://doi.org/10.1093/bioinformatics/btq228>

## Acknowledgement

I thank [Stephen Smith](https://github.com/blackrim), [Joseph Brown](https://github.com/josephwb), and [Caroline Parins-Fukuchi](https://github.com/carolinetomo) for discussions.

## License
GNU General Public License, version 3