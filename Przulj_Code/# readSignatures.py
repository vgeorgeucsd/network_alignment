# readSignatures
# Read the signatures from ndump2 files
def readSignatures(file):
    signDict = []
    
    fRead = open(file, 'r')
    
    for line in fRead:
        splitted = line.strip().split(' ')
        signDict.append([int(value) for value in splitted])
    fRead.close()
    
    return signDict