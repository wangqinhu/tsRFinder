tsRFinder
=========

A tool for tsRNA analysis with NGS data

Install
-------
```bash
git clone https://github.com/wangqinhu/tsRFinder.git
cd tsRFinder
echo export tsR_dir=$PWD >> ~/.bashrc
```

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
    -n  Min read length            [defalut 18]
    -x  Max read length            [default 45]
    -e  Min expression level       [default 10]
    -u  Mature tsRNA level cut-off [default 10]
    -f  Small RNA family threshold [default 72]
    -w  tRNA with/without label    [defualt no/yes]
    -o  Output compressed tarball  [default no/yes]
    -m  Mode, run/debug            [defualt run/debug]
    -h  Help
    -v  Version

Example:

    tsRFinder.pl -c demo/tsR.conf
```

Manual
------
See [doc/manual.pdf][1]

Demo
----
![animated gif demo][2]

[1]: doc/manual.pdf
[2]: doc/demo.gif
