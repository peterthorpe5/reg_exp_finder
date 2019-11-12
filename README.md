# reg_exp_finder
script to find reg express in biological sequences

requires biopython
pip install biopython

$ python reg_ex_finder.py.py -h
python reg_ex_finder.py.py -i in.fasta --reg expression -o outdata.txt



Options:
  -h, --help            show this help message and exit
  -i INFASTA            this is the fasta file with your sequences
  -r REG, --reg=REG     reg ex of interests
  -o FILE, --output=FILE
                        Output filename
