import json
import re
from pathlib import Path

BaseDir = Path(__file__).parent.resolve()
DataLst = ['D-RW-71.dat', 'D-RW-121.dat', 'D-RW-221.dat', 'D-RW-621.dat', 'D-RW-EXT1021.dat']

def readDRWFile(fileName):
    salt = int(re.search(r'\d+', str(fileName))[0])
    curr_seq1 = ""
    curr_seq2 = ""
    res = []
    with fileName.open() as f:
        for line in f.readlines():
            lst = line.strip().split('\t')
            if not lst[0].startswith('(re'):
                seq = lst[0].split('/')
                curr_seq1 = seq[0]
                curr_seq2 = seq[1][-1::-1]
            res.append({'seq1': curr_seq1, 'seq2': curr_seq2, 'complementary': curr_seq1 == curr_seq2, 'na_concentration': salt, 'sodium_concentration': salt-21, 'rna_concentration': float(lst[1]), 'T_exp': float(lst[2]), 'T_FIF': float(lst[3]), 'T_VIF': float(lst[4])})
    return res

def getData():
    dataPath = Path(BaseDir, 'data.json')
    if dataPath.exists():
        return json.load(dataPath.open())
    res = []
    for fileName in DataLst:
        res += readDRWFile(Path(BaseDir, fileName))

    json.dump(res, dataPath.open(mode='w'))
    return res

if __name__ == "__main__":
    getData()

