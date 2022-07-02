sp = ['Brachy','Maize','Rice','Sorghum','Telongatum']
finalBlocks_path = './example/finalBlocks'


def processGenenumber(sp, resultDir):
    block_len = {}
    for i in sp:
        with open(resultDir + '/' + i + '.final.synteny', 'r') as sf:
            for line in sf:
                temp = line
                temp = temp.rstrip('\n').rstrip()
                block = temp.split(' ')[0].split(':')[0]
                length = len(temp.split(' ')[1:])
                # print(temp.split(' ')[1:])

                if block not in block_len.keys():
                    block_len[block] = length
                else:
                    if block_len[block] < length:
                        block_len[block] = length

    with open(resultDir + '/blockindex.genenumber','w') as f:
        f.write('blockID\tblockLength\n')
        for i,j in block_len.items():
            f.write(i + '\t' + str(j) + '\n')
    return


processGenenumber(sp, finalBlocks_path)
