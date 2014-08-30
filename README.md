About
-----
The tsRFinder is a lightweight, fast and reliable tool for prediction and annotation of tRNA-derived small RNAs using next-generation sequencing data.


Installation
------------

The tsRFinder depends on:

-   Perl, v5.10.1 or higher.

-   R, v2.15.2 or higher.

-   bowtie, v1.0.0 or higher.

-   tRNAscan-SE, v1.3.1 or higher, optional.

You may clone tsRFinder by typing the following in the terminal:

    git clone https://github.com/wangqinhu/tsRFinder.git

Alternatively, you may download it from:

    https://github.com/wangqinhu/tsRFinder/releases/latest

Finally, please add tsR_dir as an environment variable:

    echo export tsR_dir="/the/path/of/tsRFinder" >> $HOME/.bashrc
    source ~/.bashrc


Usage
-----

```
tsRFinder usage:

    tsRFinder.pl <option>

    -c  Configuration file
    -l  Label
    -g  Reference genomic sequence file, conflict with -t
    -t  Reference tRNA sequence file, conflict with -g
    -s  Small RNA sequence file
    -a  Adaptor sequence
    -n  Min read length            [default 18]
    -x  Max read length            [default 45]
    -e  Min expression level       [default 10]
    -u  Mature tsRNA level cut-off [default 10]
    -f  Small RNA family threshold [default 72]
    -w  tRNA with/without label    [default no/yes]
    -o  Output compressed tarball  [default no/yes]
    -m  Mode, run/debug            [default run/debug]
    -h  Help
    -v  Version

Example:

    tsRFinder.pl -c demo/tsR.conf
```


Manual
------
A full manual could be found in [doc/manual.pdf][1]


Demo
----
![animated gif demo][2]

[1]: doc/manual.pdf
[2]: doc/demo.gif
