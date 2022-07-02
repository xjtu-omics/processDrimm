import os
from pathlib import Path
from utils import processLCSAndFirstFilter as plff
from utils import processFinalFilter as pff
import shutil


block_file = './example/drimm/blocks.txt'
drimmSyntenyFile = './example/drimm/synteny.txt'
outdir = './example'
chr_number = [5,10,12,10,7]
sp_list = ['Brachy','Maize','Rice','Sorghum','Telongatum']
target_rate = '2:4:2:2:2'

# outdir为processOrthofind的输出路径

def readSequence(file):
    sequence = []
    with open(file,'r') as f:
        while True:
            line = f.readline()[:-2]
            if not line:
                break
            itemset = line.split(' ')
            sequence.append(itemset)
    return sequence

def syntenyDict(file):
    syntenyDict = {}
    with open(file) as sf:
        while True:
            line = sf.readline()[:-2]
            if not line:
                break
            itemset = line.split(' ')
            header = itemset[0].split(':')
            syntenyDict[header[0]] = itemset[1:]
    return syntenyDict

# 用来处理drimm输出得到各个物种的输入

drimm_split_blocks_dir = outdir + '/drimmBlocks'
raw_block_dir = outdir + '/tmp'
result_dir = outdir + '/finalBlocks'

if (not Path(drimm_split_blocks_dir).exists()):
    os.makedirs(drimm_split_blocks_dir)

if (not Path(raw_block_dir).exists()):
    os.makedirs(raw_block_dir)

if (not Path(result_dir).exists()):
    os.makedirs(result_dir)


sequence = readSequence(block_file)
sp_sequences = []
last = 0
for i in range(len(chr_number)):
    sp_sequences.append(sequence[last:last+chr_number[i]])
    last += chr_number[i]

for i in range(len(sp_sequences)):
    outfile = drimm_split_blocks_dir + '/' + sp_list[i] + '.block'
    outfile = open(outfile,'w')
    for j in sp_sequences[i]:
        outfile.write('s ')
        for k in j:
            outfile.write(k+' ')
        outfile.write('\n')
    outfile.close()

processLCSAndFirstFilter = plff.processLCSAndFirstFilter(drimm_split_blocks_dir, raw_block_dir, target_rate,
                                                         drimm_split_blocks_dir, outdir, drimmSyntenyFile,
                                                         sp_list, 's')
processLCSAndFirstFilter.excute()


processFinalFilter = pff.processFinalFilter(sp_list, raw_block_dir, drimm_split_blocks_dir,  result_dir, 's')
processFinalFilter.excute()

shutil.rmtree(raw_block_dir)






# block_rate_dir = {}
# for i in sp_sequences:
#     for j in i:
#         for k in j:
#             block = ''
#             if k.startswith('-'):
#                 block = k[1:]
#             else:
#                 block = k
#             if block not in block_rate_dir.keys():
#                 rate_list = []
#                 for l in chr_number:
#                     rate_list.append(0)
#                 block_rate_dir[block] = rate_list
#
# for i in range(len(sp_sequences)):
#     for j in sp_sequences[i]:
#         for k in j:
#             block = ''
#             if k.startswith('-'):
#                 block = k[1:]
#             else:
#                 block = k
#             block_rate_dir[block][i] += 1
# save_block = []
# for i in block_rate_dir.keys():
#     rate = ''
#     for j in block_rate_dir[i]:
#         rate += str(j) + ':'
#     rate = rate[:-1]
#     if rate == target_rate:
#         save_block.append(i)
#
# synteny = syntenyDict(synteny_file)
# save_block_filter = []
# for i in save_block:
#     save_block_filter.append(i)
# # 输出过滤情况
# print(len(save_block))
# print(len(save_block_filter))
#


