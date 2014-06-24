skewer
======

<b>skewer</b> (transferred from https://sourceforge.net/projects/skewer) implements the <i>bit-masked k-difference matching algorithm</i> dedicated to the task of adapter trimming and it is specially designed for processing next-generation sequencing (NGS) paired-end sequences.

### Citation
Jiang, H., Lei, R., Ding, S.W. and Zhu, S. (2014) Skewer: a fast and accurate adapter trimmer for next-generation sequencing paired-end reads. <i>BMC Bioinformatics</i>, <b>15</b>, 182.

### Features
* Detection and removal of adapter sequences
* Insertion and deletion allowed in pattern matching
* Targeted at Single End, Paired End (PE), and Long Mate Pair (LMP) reads
* Demultiplexing of barcoded sequencing runs
* Multi-threading support
* Trimming based on phred quality scores
* IUPAC characters for barcodes and adapters
* Compressed input and output support

### Installation from binary
Copy <b>skewer</b> to your favorate BIN directory, and make sure the PATH environment variable is correctly set. For example:

    $ mkdir -p ~/bin
    $ cp -p skewer ~/bin/
    $ echo 'export PATH=~/bin:$PATH' >> ~/.bashrc
    $ source ~/.bashrc

### Installation from source codes
Enter into the directory of source codes, then:

    $ make
    $ sudo make install
