# MTPT全步骤分析程序

## 用到的库

Bio库和pandas库


## 参数列表

- cp 叶绿体基因的GenBank文件
- mt 线粒体的GenBank文件
- blastnfile balstn比对结果txt文件，利用***线粒体***建库，和叶绿体进行比对
- length 筛选的最小长度
- scores 筛选的最低得分
>例如 `python MTPT.py  --cp .\Chloranthus_spicatus_cp.gb --mt .\Chloranthus_spicatus_mt.gb --blastnfile .\Chloranthus_spicatus_result.txt --length 50 --scores 80`

### 如何生成blastnfile文件
1. 利用线粒体建库
- makeblastdb -in ***线粒体基因文件*** -dbtype nucl -parse_seqids -out ***线粒体数据库***
2. 叶绿体和线粒体比对
- blastn -query  ***叶绿体基因文件*** -db ***线粒体数据库***  -out ***结果文件*** -task blastn -evalue 1e-5 -word_size 9 -gapextend 2 -reward 2 -penalty -3 -gapopen 5 -outfmt 6 -num_threads 4

## 生成结果说明

 会先生成一个result文件夹，里面会包含以下文件：
- MTPT.xlsx MTPT分析结果表格
- all_seq_length.txt 线粒体和叶绿体的整个序列长度
- mtpt_seq_length.txt 转移序列长度
- GC_content.txt GC含量
- escape_seq.txt 转移序列

  
