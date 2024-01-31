# _*_ coding :utf-8 _*_
# @Time :2024/1/29
# @Author: he
# @File: MTPT
# @Project: MTPT

import os
from Bio import SeqIO
import sys
import getopt
import pandas as pd

result_path = 'result/'


# 计算GC含量
def get_GC_content(write, seq, name):
    for i in range(0, len(seq), 1000):
        sub_seq = seq[i:i + 1000]
        start = i
        end = i + len(sub_seq)
        gc_num = sub_seq.count('G') + sub_seq.count('C')
        gc = gc_num / len(sub_seq)
        gc = round(gc, 2)
        write.write(name + "\t" + str(start) + "\t" + str(end) + "\t" + str(gc) + "\n")


# 生成GC含量文件
def set_GC_content(cp_file, mt_file):
    cp = SeqIO.read(cp_file, 'genbank')
    mt = SeqIO.read(mt_file, 'genbank')
    cp_seq = cp.seq
    mt_seq = mt.seq
    file = result_path + 'GC_content.txt'
    print(file)
    result = open(file, 'w')
    get_GC_content(result, cp_seq, 'cpDNA')
    get_GC_content(result, mt_seq, 'mtDNA')
    result.close()


# 生成序列总长度文件
def set_all_seq_length(cp_file, mt_file):
    cp = SeqIO.read(cp_file, 'genbank')
    mt = SeqIO.read(mt_file, 'genbank')
    cp_seq = cp.seq
    mt_seq = mt.seq
    file = result_path + 'all_seq_length.txt'
    result = open(file, 'w')
    result.write('cpDNA\t' + str(len(cp_seq)) + '\t251,180,174\n')
    result.write('mtDNA\t' + str(len(mt_seq)) + '\t143,178,200\n')


# 得到比对结果
def get_blastn(blastn_file, min_length, min_scores):
    blastns = []
    with open(blastn_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            identity = line[2]
            alignment_length = line[3]
            mismatches = line[4]
            gap_openings = line[5]
            cp_start = line[6]
            cp_end = line[7]
            mt_start = line[8]
            mt_end = line[9]
            blsatn = {
                'identity': eval(identity),
                'alignment_length': eval(alignment_length),
                'mismatches': eval(mismatches),
                'gap_openings': eval(gap_openings),
                'cp_start': eval(cp_start),
                'cp_end': eval(cp_end),
                'mt_start': eval(mt_start),
                'mt_end': eval(mt_end)
            }
            if eval(alignment_length) >= min_length and eval(identity) >= min_scores:
                blastns.append(blsatn)
    return blastns


# 生成转移序列长度文件
def set_seq_length(blastns):
    seq_result = open(result_path + 'mtpt_seq_length.txt', 'w')
    for blsatn in blastns:
        cp_start = blsatn['cp_start']
        cp_end = blsatn['cp_end']
        mt_start = blsatn['mt_start']
        mt_end = blsatn['mt_end']
        seq_length = blsatn['alignment_length']
        seq_result.write(f'cpDNA\t{cp_start}\t{cp_end}\t{seq_length}\n')
        seq_result.write(f'mtDNA\t{mt_start}\t{mt_end}\t{seq_length}\n')
    seq_result.close()


# 生成转义序列线图文件
def set_escape_seq(blastns):
    escape = open(result_path + 'escape_seq.txt', 'w')
    for blsatn in blastns:
        cp_start = blsatn['cp_start']
        cp_end = blsatn['cp_end']
        mt_start = blsatn['mt_start']
        mt_end = blsatn['mt_end']
        mtpt_length = blsatn['alignment_length']
        if mtpt_length >= 1000:
            escape.write(f'cpDNA\t{cp_start}\t{cp_end}\tmtDNA\t{mt_start}\t{mt_end}\t38,70,83\n')
        else:
            escape.write(f'cpDNA\t{cp_start}\t{cp_end}\tmtDNA\t{mt_start}\t{mt_end}\t42,157,142\n')
    escape.close()


# 生成execl表格
def set_excel(cp_file, mt_file, blastn_file):
    data = {
        'No.': [],
        'Alignment length (bp)': [],
        'Identity (%)': [],
        'Mis-matchs': [],
        'Gap openings': [],
        'CP start': [],
        'CP end': [],
        'MT start': [],
        'MT end': [],
        'MTPT annotaton': []
    }
    cp_gene = parse(cp_file)
    mt_gene = parse(mt_file)
    index = 1
    for line in blastn_file:
        cp_start = line['cp_start']
        cp_end = line['cp_end']
        mt_start = line['mt_start']
        mt_end = line['mt_end']
        cp_result, mt_result = get_gene_in_mtpt(line, cp_gene, mt_gene)
        partial = set()
        complete = set()
        for gene in cp_result:
            if 'partial' in gene:
                partial.add(gene.split('_')[0])
            if 'complete' in gene:
                complete.add(gene.split('_')[0])
        for gene in mt_result:
            if 'partial' in gene:
                partial.add(gene.split('_')[0])
            if 'complete' in gene:
                complete.add(gene.split('_')[0])
        partial_str = ''
        complete_str = ''
        if len(partial) != 0:
            partial_str = 'Partial(' + ','.join(partial) + ') '
        if len(complete) != 0:
            complete_str = 'Complete(' + ','.join(complete) + ')'
        data['No.'].append('MTPT' + str(index))
        index += 1
        data['Alignment length (bp)'].append(line['alignment_length'])
        data['Identity (%)'].append(line['identity'])
        data['Mis-matchs'].append(line['mismatches'])
        data['Gap openings'].append(line['gap_openings'])
        data['CP start'].append(cp_start)
        data['CP end'].append(cp_end)
        data['MT start'].append(mt_start)
        data['MT end'].append(mt_end)
        data['MTPT annotaton'].append(partial_str + complete_str)
    df = pd.DataFrame(data)
    # 将DataFrame写入Excel文件
    df.to_excel(result_path + 'MTPT.xlsx', index=False)


# 得到命令行参数的文件路径
def get_args():
    opts, args = getopt.getopt(sys.argv[1:], '-h-v-c-m-b-l-s:v',
                               ['help', 'cp=', 'mt=', 'version', 'blastnfile=', 'length=', 'scores='])
    min_length, min_scores = 0, 0
    for opt_name, opt_value in opts:
        if opt_name in ('-h', '--help'):
            print("[*] Help info")
            print("[-] Usage: python MTPT.py --cp <叶绿体基因的GenBank文件> --mt <线粒体的GenBank文件> --blastnfile "
                  "<利用叶绿体建库的叶绿体和线粒体的balstn比对结果txt文件> --length <最小长度> --scores <最低得分>")
            sys.exit(0)
        if opt_name in ('-v', '--version'):
            print("[*] Version is 0.01 ")
            sys.exit(0)
        if opt_name in ('-c', '--cp'):
            cp = opt_value
        if opt_name in ('-m', '--mt'):
            mt = opt_value
        if opt_name in ('-b', '--blastnfile'):
            blastnfile = opt_value
        if opt_name in ('-l', '--length'):
            min_length = eval(opt_value)
        if opt_name in ('-s', '--scores'):
            min_scores = eval(opt_value)
    return cp, mt, blastnfile, min_length, min_scores


# 解析GenBank文件
def parse(genebank):
    gb_result = SeqIO.parse(open(genebank, errors='ignore'), 'genbank')
    record = next(gb_result, None)
    features = record.features
    gb = []
    for f in features:
        if f.type == 'CDS':
            try:
                gene_name = f.qualifiers['gene'][0]
            except KeyError:
                continue
            location = f.location
            gb.append({'gene_name': gene_name, 'location': location})
        if f.type == 'rRNA':
            try:
                gene_name = f.qualifiers['gene'][0]
            except KeyError:
                continue
            location = f.location
            gb.append({'gene_name': gene_name, 'location': location})
        if f.type == 'tRNA':
            try:
                gene_name = f.qualifiers['gene'][0]
            except KeyError:
                continue
            location = f.location
            gb.append({'gene_name': gene_name, 'location': location})
    return gb


# 得到在转移序列上的基因
def get_gene_in_mtpt(mtpt_line, cp_gene, mt_gene):
    cp_start = mtpt_line['cp_start']
    cp_end = mtpt_line['cp_end']
    cp_result = []
    for gene in cp_gene:
        gene_name = gene['gene_name']
        location = gene['location']
        parts = location.parts
        gene_start = parts[0].start
        gene_end = parts[-1].end
        if cp_start <= gene_start <= gene_end <= cp_end:
            cp_result.append(gene_name + '_complete')
        if gene_start <= cp_start <= gene_end <= cp_end or cp_start <= gene_start <= cp_end <= gene_end:
            cp_result.append(gene_name + '_partial')
    mt_start = mtpt_line['mt_start']
    mt_end = mtpt_line['mt_end']
    mt_result = []
    for gene in mt_gene:
        gene_name = gene['gene_name']
        location = gene['location']
        parts = location.parts
        gene_start = parts[0].start
        gene_end = parts[-1].end
        if mt_start <= gene_start <= gene_end <= mt_end:
            mt_result.append(gene_name + '_complete')
        if gene_start <= cp_start <= gene_end <= cp_end or cp_start <= gene_start <= cp_end <= gene_end:
            mt_result.append(gene_name + '_partial')
    return cp_result, mt_result


if __name__ == '__main__':
    cp_file, mt_file, blastn_file, length, scores = get_args()
    if not os.path.exists(result_path):
        os.makedirs(result_path)
    set_GC_content(cp_file, mt_file)
    set_all_seq_length(cp_file, mt_file)
    blastn_result = get_blastn(blastn_file, length, scores)
    set_seq_length(blastn_result)
    set_escape_seq(blastn_result)
    set_excel(cp_file, mt_file, blastn_result)
