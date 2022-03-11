#!/bin/sh

cd asap2/

# taxonomic classifier weighted Silva 138 
python -u ../bin/asap2.py -i data -c silva-138-99-nb-weighted-classifier.qza -p 10 -o maleth
