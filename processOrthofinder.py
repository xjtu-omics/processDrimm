import pandas as pd
import numpy as np

dir = './example'
sp = ['Brachy','Maize','Rice','Sorghum','Telongatum']
gff_list = ['Brachy.gff','Maize.gff','Rice.gff','Sorghum.gff','Telongatum.gff']
sp_ratio = [2,4,2,2,2]



def getFilterSequence(bed, group_filter_dir):
    bed = pd.read_csv(bed,sep='\t',header=None)[[0,1]]

    chrlist = bed[0].unique().tolist()

    new_sequence = []
    for i in chrlist:
        split = bed.loc[bed[0] == i][1].tolist()

        chr = []
        for j in split:
            if j in group_filter_dir.keys():
               chr.append(group_filter_dir[j])
        new_sequence.append(chr)
    return new_sequence

def getAllSequence(bed, group_dir):
    bed = pd.read_csv(bed,sep='\t',header=None)[[0,1]]

    chrlist = bed[0].unique().tolist()

    new_sequence = []
    new_name_sequence = []
    for i in chrlist:
        split = bed.loc[bed[0] == i][1].tolist()

        chr = []
        chr_name = []
        for j in split:
            if j in group_dir.keys():
               chr.append(group_dir[j])
               chr_name.append(j)

        new_sequence.append(chr)
        new_name_sequence.append(chr_name)
    return new_sequence,new_name_sequence


def outSequence(sequence,outfile):
    outfile = open(outfile,'w')
    for i in sequence:
        for j in i:
            outfile.write(str(j)+' ')
        outfile.write('\n')
    outfile.close()


gff_path_list = []
for i in range(len(gff_list)):
    temp = dir + '/' + gff_list[i]
    gff_path_list.append(temp)

for i in gff_path_list:
    dataFrame = pd.read_csv(i, header=None, sep='\t')
    name = dataFrame.at[0, 0]
    number_index = 0
    for j in range(len(name) - 1, -1, -1):
        if name[j].isdigit():
            continue
        else:
            number_index = j + 1
            break
    # print(number_index)
    # print(name[number_index:])
    dataFrame[0] = dataFrame[0].map(lambda x: int(x[number_index:]))

    dataFrame = dataFrame.sort_values(by=[0, 2], ascending=[True, True])
    dataFrame[0] = dataFrame[0].map(lambda x: name[:number_index] + str(x))
    dataFrame[1] = dataFrame[1].map(lambda x: x.split(';')[0])
    dataFrame.to_csv(i, header=None, index=None, sep='\t')

ortho = pd.read_csv(dir+'/Orthogroups.tsv',sep='\t')
ortho = ortho.fillna('')
# print(ortho)
columns = ortho.columns.tolist()


ortho = np.asarray(ortho)
# print(ortho)
group_dir = {}
for i in ortho:
    group = i[0]
    group_dir[group] = {}
    species = i[1:]
    for j in range(len(species)):
        genes = species[j].split(', ')
        if genes[0] == '':
            group_dir[group][sp[j]] = []
        else:
            group_dir[group][sp[j]] = genes

rate_dir = {}
finalGroup = {}

for i in group_dir.keys():
    rate_list = []
    for j in group_dir[i].keys():
        rate_list.append(len(group_dir[i][j]))
    ok = 1
    for j in range(len(rate_list)):
        if rate_list[j] > sp_ratio[j] or rate_list[j] == 0:
            ok = 0
    if ok == 0:
        continue
    else:
        rate = ''
        for j in rate_list:
            rate += str(j) + ':'
        rate = rate[:-1]
        finalGroup[i] = group_dir[i]
        if rate not in rate_dir.keys():
            rate_dir[rate] = 1
        else:
            rate_dir[rate] += 1
print('gene rate')
for i in rate_dir.keys():
    print(i + '\t' + str(rate_dir[i]))

outfile = dir + '/group.xls'
count = 1
outfile = open(outfile,'w')
outfile.write('gene\tgroup\n')

outfile_filter = dir + '/filter_group.xls'
outfile_filter = open(outfile_filter,'w')
outfile_filter.write('gene\tgroup\n')

for i in group_dir.keys():
    for j in group_dir[i].keys():
        for k in group_dir[i][j]:
            outfile.write(k+'\t'+str(count)+'\n')
    if i in finalGroup.keys():
        for j in finalGroup[i].keys():
            for k in finalGroup[i][j]:
                outfile_filter.write(k + '\t' + str(count) + '\n')
    count += 1
outfile.close()
outfile_filter.close()

group = pd.read_csv(dir +'/group.xls',sep='\t')
group = np.asarray(group)
group_dir = {}
for i in group:
    group_dir[i[0]] = i[1]

group_filter = pd.read_csv(dir +'/filter_group.xls',sep='\t')
group_filter = np.asarray(group_filter)
group_filter_dir = {}
for i in group_filter:
    group_filter_dir[i[0]] = i[1]

sample_sequence_files = dir + '/drimm.sequence'
sample_sequence_files = open(sample_sequence_files, 'w')

for i in gff_list:
    gff = dir + '/' + i
    sequence,sequence_name = getAllSequence(gff, group_dir)
    filter_sequence = getFilterSequence(gff, group_filter_dir)
    item = i.split('.')
    outfile = dir  + '/' + item[0] +'.sequence'
    outSequence(filter_sequence, outfile)
    outallfile = dir + '/' + item[0] +'.all.sequence'
    outSequence(sequence, outallfile)
    outallfilename = dir + '/' + item[0] + '.all.sequence.genename'
    outSequence(sequence_name, outallfilename)

    for j in filter_sequence:
        for k in j:
            sample_sequence_files.write(str(k) + ' ')
        sample_sequence_files.write('\n')

sample_sequence_files.close()
