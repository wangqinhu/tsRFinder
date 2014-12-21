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

If you want to run tsRFinder as a system command, create a soft link for it:

	ln -s `pwd`/tsRFinder/tsRFinder.pl /usr/local/bin/

Usage
-----

A typical tsRFinder job can be easily finished on a modern MacBook (running OS X) or laptop (running Linux), see the following usage and examples.

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
    -i  interactive                [default yes/no]
    -m  Mode, run/debug            [default run/debug]
    -h  Help
    -v  Version

Examples:

    tsRFinder.pl -c demo/tsR.conf

    tsRFinder.pl -c demo/tsR.alt.conf
```

To submit tsRFinder job on a cluster managed with sun grid engine (sge), see [demo/tsR.sge.sh][3]

Example:

	qsub demo/tsR.sge.sh demo/tsR.conf


Manual
------
A full manual could be found in [doc/manual.pdf][1]


Demo
----
![animated gif demo][2]

[1]: https://raw.githubusercontent.com/wangqinhu/tsRFinder/master/doc/manual.pdf
[2]: https://raw.githubusercontent.com/wangqinhu/tsRFinder/master/doc/demo.gif
[3]: https://raw.githubusercontent.com/wangqinhu/tsRFinder/master/demo/tsR.sge.sh
