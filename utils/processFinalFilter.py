from copy import deepcopy


class processFinalFilter:
    def __init__(self, sp, blockDir, syntenyDir, resultDir, chr_shape='s'):
        self.sp = sp
        self.blockDir = blockDir
        self.syntenyDir = syntenyDir
        self.resultDir = resultDir
        self.chr_shape = chr_shape

    def outputFilterBlock(self, blockDir, syntenyDir, resultDir, sp):

        # 获取synteny所有匹配成功的block以及次数
        synteny_block_position = {}
        with open(syntenyDir + '/' + sp + '.synteny', 'r') as sf:
            for i in sf:
                block_info = i
                block_info = block_info.split(' ')[0].split(':')
                if block_info[0] not in synteny_block_position.keys():
                    synteny_block_position[block_info[0]] = []
                synteny_block_position[block_info[0]].append(int(block_info[1]))

        # 获取所有block序列
        block_sequence = []
        with open(blockDir + '/' + sp + '.block', 'r') as inf:
            for i in inf:
                block_info = i
                block_info = block_info.rstrip('\n').rstrip()
                block_info = block_info.split(' ')[1:]
                block_sequence.append(block_info)

        # 初始化block出现序列为0
        raw_block_position = {}
        for i in block_sequence:
            for j in i:
                if j.startswith('-'):
                    j = j[1:]
                if j not in raw_block_position.keys():
                    raw_block_position[j] = 0
        block_list = deepcopy(raw_block_position)

        with open(resultDir + '/' + sp + '.final.block', 'w') as outf:
            for i in block_sequence:
                if self.chr_shape == 's' or self.chr_shape == 'S':
                    outf.write('s ')
                else:
                    outf.write('c ')

                for j in i:
                    if j.startswith('-'):
                        block = j[1:]
                    else:
                        block = j

                    raw_block_position[block] += 1
                    if raw_block_position[block] in synteny_block_position[block]:
                        outf.write(j + ' ')
                    # else:
                    #     print(block+':'+str(raw_block_position[block]) + 'drop!')
                outf.write('\n')
        return block_list

    def outputFilterSynteny(self, block_list, syntenyDir, resultDir, sp):
        block_list_copy = deepcopy(block_list)
        with open(resultDir + '/' + sp + '.final.synteny', 'w') as synWriteFile:
            with open(syntenyDir + '/' + sp + '.synteny', 'r') as synReadFile:
                for i in synReadFile:
                    temp = i.rstrip('\n').rstrip()
                    temp = temp.split(' ')
                    block = temp[0].split(':')[0]
                    if block in block_list.keys():
                        block_list[block] += 1
                        synWriteFile.write(block + ':' + str(block_list[block]) + ':' + temp[0].split(':')[2] + ':'  + temp[0].split(':')[3] + ' ')
                        synWriteFile.write(' '.join(temp[1:]) + ' \n')

        with open(resultDir + '/' + sp + '.final.synteny.genename', 'w') as synWriteFile:
            with open(syntenyDir + '/' + sp + '.synteny.genename', 'r') as synReadFile:
                for i in synReadFile:
                    temp = i.rstrip('\n').rstrip()
                    temp = temp.split(' ')
                    block = temp[0].split(':')[0]
                    if block in block_list_copy.keys():
                        block_list_copy[block] += 1
                        synWriteFile.write(block + ':' + str(block_list_copy[block]) + ':'  + temp[0].split(':')[2] + ':'  + temp[0].split(':')[3] + ' ')
                        synWriteFile.write(' '.join(temp[1:]) + ' \n')

        return

    def excute(self):

        for i in self.sp:
            temp = self.outputFilterBlock(self.blockDir, self.syntenyDir, self.resultDir, i)
            self.outputFilterSynteny(temp, self.syntenyDir, self.resultDir, i)

# if __name__ == '__main__':
#     sp = ['Brachy', 'Maize', 'Rice', 'Sorghum','Telongatum']
#     c = processFinalFilter(sp,'C:/Users/guo/OneDrive/IAGS_input_builder/IAGS_input_builder/outputFinal/temp/firstFilterBlock',
#                            'C:/Users/guo/OneDrive/IAGS_input_builder/IAGS_input_builder/outputFinal/result/drimmBlock',
#                            'C:/Users/guo/OneDrive/IAGS_input_builder/IAGS_input_builder/outputFinal/result/finalBlock','s')
#     c.excute()
