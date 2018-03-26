# Arachis

Arachis is a python 3.* package that allow users to reconstruct ancestral genome gene orders and 
infer pairwise plastome differences or events.

The algorithm for reconstructing ancestral plastome gene orders implement in script file `run_pypmag.py`
is derived from the ancestral gene order reconstruction module of `PMAG+`, with modifications:
1. Accept circular and gap-containing genomes as inputs. See modification on GRIMM format below.
2. Equipped with python multiprocessing.

This package implement a derived version of [the classic GRIMM format](http://grimm.ucsd.edu/GRIMM/grimm_instr.html) 
with following modifications:
1. Blocks could be named with letters in "`-.0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ_abcdefghijklmnopqrstuvwxyz|~[]`". 
"`-`" still means a reverse direction when it appears at the first letter.
2. The "`*`" block in a sequence means a gap. `1 2 3 * 4 5 * 6 $` is equivalent to `1 2 3 * -5 -4 * 6 $`.
3. A sequence line without the "`$`" block in the end stands for a circular chromosome. But in this case only 
one chromosome per sample is allowed, and multiple lines without "$" would be regarded as one single chromosome 
written interleaved. This design is due to the limitation of applying to tsp solver.

The functions in Arachis for inferring pairwise plastome differences or events are still at <b>infant</b> stage. 
If you find any bugs or something to improve, please contact 
[jinjianjun@mail.kib.ac.cn](mailto:jinjianjun@mail.kib.ac.cn). New contributors are welcome!

Also, users have to bear in mind that do not test data with too many breakpoints (like 10+).
The function SignedPermutation.inversion_event_from utilizes an <b>exhausted</b> 
scheme searching for one best solution.

## Installation

Download Arachis and install Arachis with:

    git clone ...
    cd ARACHIS
    python setup.py install

To use `run_pypmag.py` to reconstruct ancestral genome gene order, you have to install following dependencies:

1. RAxML
2. tsp_solver

## Example
Run `pypmag` with test data:

    python run_pypmag.py -d test/test_1_grimm.txt -t test_1_grimm.tre -o test/test_1_output --seed 12345
    
To see parsimonious events along the branch from `A1` to `sp2`:

    python
    >>> from arachis.permutationClass import GenomeList
    >>> extant_samples = GenomeList("test/test_1_grimm.txt")
    >>> sp2 = extant_samples["sp2"].chromosomes()[0]
    >>> ancestors = GenomeList("test/test_1_output/OutputOrder")
    >>> A1 = ancestors["A1"].chromosomes()[0]
    >>> events = sp2.event_from(A1)

## Citation
If you use Arachis in your research, please cite this github page.
If you use run_pypmag.py, please cite following papers as well:
1. RAxML
2. tsp_solver
3. pmag+
4. dendropy

## License
GNU General Public License, version 3