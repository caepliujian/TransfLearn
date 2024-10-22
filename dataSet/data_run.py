from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys, AllChem
from rdkit.Chem.PandasTools import pd
from tqdm import tqdm


def getTestSample(inputFile,n):
    data = pd.read_csv(inputFile, delim_whitespace=False, header=0)
    test = data.sample(n=n)
    # test = test.drop(labels=test.index)
    test.to_csv(r"testSet.csv")
    # data = data.drop(labels=test.index)
    # data.to_csv(r"./data/ttt.csv")


def mvRow(file1, file2, out):
    lines1 = open(file1, "r").readlines()
    lines2 = open(file2, "r").readlines()
    lines1 = set(lines1)
    lines2.pop(0)
    lines2 = set(lines2)
    result = lines1 - lines2
    open(out, "x").writelines(result)
    print(len(lines1), len(lines2), len(result))


def dataStatistics(file, pc_line, d_line, out):
    lines = open(file, "r").readlines()
    header = lines.pop(0)
    writer = open(out, "x")
    writer.write(header)
    counter = 0
    result = []
    for line in tqdm(lines):
        values = line.split(",")
        density = float(values[2])
        smiles = Chem.CanonSmiles(values[1])
        pc = float(values[5])
        if density > d_line or pc > pc_line:
            writer.write(line)
            counter += 1
            result.append(line)
    print(counter)
    return result


def getDifferentMol(allMolFile, trainFile, outFile):
    testSet = open(r"C:\Users\tangyuechuan\Desktop\chemProp\data\testSet\testSet.csv", "r").readlines()
    testSet = set(testSet)
    allList = open(allMolFile, "r").readlines()
    trainList = pd.read_csv(trainFile, delim_whitespace=False, header=0)
    trainSmilesList = trainList["smiles"]
    # AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), 2)
    # MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smiles))
    # trainMolFP = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), 2) for smiles in trainSmilesList]
    trainMolFP = [MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smiles)) for smiles in trainSmilesList]
    header = allList.pop(0)
    smilesFind = []
    for line in tqdm(allList):
        smiles = line.split(",")[1]
        # baseFP = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), 2)
        baseFP = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smiles))
        isFine = True
        for fp in trainMolFP:
            similar = DataStructs.FingerprintSimilarity(fp, baseFP)
            if similar > 0.7:
                isFine = False
                break
        if isFine:
            if line not in testSet:
                trainMolFP.append(baseFP)
                smilesFind.append(line)
    writer = open(outFile, "x")
    writer.write(header)
    writer.writelines(smilesFind)


def rmExcessive(path,out):
    lines = open(path, "r").readlines()
    header = lines.pop(0)
    nameList = set()
    row = []
    for line in lines:
        name = line.split(",")[0]
        if name in nameList: continue
        row.append(line)
        nameList.add(name)
    writer = open(out,"x")
    writer.write(header)
    writer.writelines(row)
    print(len(lines), len(row)+1)


def calcSMForEach(inputFile):
    lines = open(inputFile,"r").readlines()
    lines.pop(0)
    avgSM = []
    # MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smiles))
    # AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(line.split(",")[1]), 3)
    molFPList = [MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(line.split(",")[1])) for line in lines]
    total = len(molFPList)
    for fp1 in molFPList:
        sm = 0
        for fp2 in molFPList:
            sm += DataStructs.FingerprintSimilarity(fp2,fp1)
        avgSM.append((sm-1)/total)
    for i in range(total):
        print(avgSM[i],lines[i].split(",")[1])


def getMolFromSke(ske,molFile):
    data = pd.read_csv(molFile, delim_whitespace=False, header=0)
    smilesList = data["smiles"]
    skeMol = Chem.MolFromSmarts(ske)
    for smiles in smilesList:
        mol = Chem.MolFromSmiles(smiles)
        if mol.HasSubstructMatch(skeMol):
            print(smiles)

def devCount():
    lines = open(r"C:\Users\tangyuechuan\Desktop\chemProp\predict\all_nitro_pred.csv","r").readlines()
    lines.pop(0)
    for k in range(20):
        devCount = 0
        num = 0
        n = 0.5
        start = k*0.015+n
        end = (k+1)*0.015+n
        for line in lines:
            values = line.split(",")
            d = float(values[6])
            pc = float(values[5])
            dev = float(values[7])
            if start<pc<end:
                num += 1
                devCount += dev
        print("%.4f<d<%.4f"%(start,end),num,devCount/num)
if __name__ == '__main__':
    # getTestSample(r"./orgSet/d_and_pc.csv",100)
    # mvRow(r"C:\Users\tangyuechuan\Desktop\chemProp\data\trainSet\magnet_org2002.csv"
    #       , r"C:\Users\tangyuechuan\Desktop\chemProp\data\testSet\testSet.csv"
    #       , r"C:\Users\tangyuechuan\Desktop\chemProp\data\trainSet\magnet.csv")
    # rmExcessive(r"C:\Users\tangyuechuan\Desktop\chemProp\data\trainSet\pc_d\pc_d_trainSet.csv"
    #             ,r"C:\Users\tangyuechuan\Desktop\chemProp\data\trainSet\pc_d\pc_d_trainSet——2.csv")
    # dataStatistics(r"C:\Users\tangyuechuan\Desktop\chemProp\data\orgSet\nitro_Set.csv"
    #                , 0.72, 1.7
    #                ,"./orgSet/d_or_pc.csv")
    # getDifferentMol(r"C:\Users\tangyuechuan\Desktop\chemProp\data\trainSet\nitro_ts.csv"
    #                 ,r"C:\Users\tangyuechuan\Desktop\chemProp\data\trainSet\d\d_a_non_nitro.csv"
    #                 ,r"C:\Users\tangyuechuan\Desktop\chemProp\data\trainSet\d\dis_similar_MACCS_0.7.csv")
    # calcSMForEach(r"C:\Users\tangyuechuan\Desktop\chemProp\data\orgSet\d_1.8.csv")
    # getMolFromSke("c1(nnn2)c2cnnc1",r"D:\pycharmProject\chemprop-master\molSplit\genRes\csv\genRes.csv")
    # getMolFromSke("[N+](=O)[O-]",r"C:\Users\tangyuechuan\Desktop\chemProp\data\orgSet\nitro_Set.csv")
    # devCount()
    data = pd.read_csv(r"C:\Users\tangyuechuan\Desktop\chemProp\data\test.csv", delim_whitespace=False, header=0)
    s = data["smiles"][1]
    print(s)
    print([s,r"\ff"])
    data = {'name': [s]}
    print(data)