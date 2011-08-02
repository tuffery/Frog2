import time
from frowns.parsers import smiles_parser
from frowns.Canonicalization import Traverse

from test_smiles import smilesStrings

def test():
    start = time.time()
    for smile in smilesStrings:
        print '*'*44

        tokens = smiles_parser.tokenizer(smile)
        for event, dict in tokens:
            for key, val in dict.items():
                if val == None:
                    del dict[key]

            print '\t',event, dict
    end = time.time()
    print (end-start)/len(smilesStrings), "per smile"
    print "%s seconds %s smiles"%((end-start), len(smilesStrings))

test()
