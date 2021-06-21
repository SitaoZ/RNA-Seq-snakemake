#coding:utf-8
#!/usr/bin/python

import os
import re
import sys
import json
import argparse

## python 3.6 supported
## modefied 20210528
## Auto generate analysis shell script

def parse_args():
    parser = argparse.ArgumentParser(description='RNAseq analysis pipeline')
    parser.add_argument('--samples', '-s',type=str, help='RNAseq sample list')
    parser.add_argument('--config', '-c',type=str, help='Configure path for analysys')
    parser.add_argument('--output', '-o',type=str, help='output directory')
    return parser.parse_args()

def sample_parse(sample_path):
    """
    Sample fastq file parse
    :param Sample_Path: the file contain sample's fastq file absolute path
    :return: dict :{'sample name':['fastq1','fastq2']}
    """
    samples = dict()
    with open(sample_path, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if not len(line) or line.startswith('#'):
                continue
            else:
                ss = line.split('\t')
                samples[ss[0]] = [ss[1], ss[2]]

    return samples


def get_options(text):
    matches=re.findall(r'\"(.+?)\"',text)
    print("matches:",matches)
    return ','.join(matches)

def get_option_key(text):
    text = text.replace(' ','') # remove space
    return text

def configure(configure_path):
    """
    Configure file for RNAseq analysis
    :param Configure_Path:
    :return: dict config
    """
    with open(configure_path,'r') as f:
        conf_dict = json.load(f)

    return conf_dict


def check_dir(path):
    """
    Create analysis directory
    :param path: The path for directort
    :return:
    """
    if not os.path.exists(path):
        os.makedirs(path)
        return path
    else:
        print (path,"the path exists !")
        return path






## step1 filter
def filter_cutadapt(samples, conf_dict, result_dir):
    samples_dict = dict()
    filter_shell_dir = os.path.join(result_dir, "shell/filter_cutadapt")
    filter_result_dir = os.path.join(result_dir, "result/filter_cutadapt")
    cutadapt = conf_dict['cutadapt']
    for i in samples.keys():
        singel_filter_shell_dir = check_dir(os.path.join(filter_shell_dir,i))
        singel_filter_result_dir = check_dir(os.path.join(filter_result_dir,i))

        shell_file = os.path.join(filter_shell_dir, i, "work.sh")
        with open(shell_file,'w') as shell:
            shell.writelines("zcat %s %s |gzip > %s/%s.fq.gz\n"% \
                             (samples[i][0],samples[i][1],singel_filter_shell_dir,i))
            shell.writelines("%s %s %s/%s.fq.gz\n" % \
                             (cutadapt['path'],cutadapt['options'],\
                              singel_filter_result_dir,i))
        samples_dict[i] = "%s/%s.fq.gz" % (singel_filter_result_dir,i)
        
    return samples_dict


## step2 alignment
##  hisat alignment of genome reference
def genomeMapping_hisat(conf_dict, samples_dict, result_dir):
    hisat = conf_dict['hisat']
    hisat_path = hisat['path']
    hisat_options_pe = hisat['hisat_options_pe']
    hisat_index = hisat['hisat_index']

    hisat_shell_dir = os.path.join(result_dir, "shell/genome_mapping_hisat")
    hisat_result_dir = os.path.join(result_dir, "result/genome_mapping_hisat")

    for i in samples_dict.keys():
        singel_hisat_shell_dir = check_dir(os.path.join(hisat_shell_dir, i))
        singel_hisat_result_dir = check_dir(os.path.join(hisat_result_dir, i))

        shell_file = os.path.join(singel_hisat_shell_dir, "work.sh")
        with open(shell_file, 'w') as shell:
            shell.writelines("%s %s -x %s -U %s -S %s.sam 2>%s.stat\n" % \
                             (hisat_path,hisat_options_pe,hisat_index,samples_dict[i],\
                              singel_hisat_result_dir,singel_hisat_result_dir))
            shell.writelines("%s view %s.sam -b -S -o %s.bam\n" % \
                             (conf_dict['samtools'],singel_hisat_result_dir,singel_hisat_result_dir))



# bowtie alignment of gene reference (for gene expression analysis)
def gene_mapping_bowtie(conf_dict, samples_dict, result_dir):
    bam_dict = {}
    bowtie = conf_dict['bowtie']
    bowtie_path = bowtie['path']
    genebowtie2index = bowtie['bowtie2index']
    genebowtie2_options_pe = bowtie['genebowtie2_options_pe']

    bowtie_shell_dir = os.path.join(result_dir, "shell/gene_mapping_bowtie")
    bowtie_result_dir = os.path.join(result_dir, "result/gene_mapping_bowtie")
    for i in samples_dict.keys():
        singel_bowtie_shell_dir = check_dir(os.path.join(bowtie_shell_dir, i))
        singel_bowtie_result_dir = check_dir(os.path.join(bowtie_result_dir, i))
        shell_file = os.path.join(singel_bowtie_shell_dir, "work.sh")
        with open(shell_file,'w') as shell:
            shell.writelines("%s %s -x %s -U %s > %s/%s.sam\n" % \
                             (bowtie_path, genebowtie2_options_pe, \
                              genebowtie2index, samples_dict[i],\
                              singel_bowtie_result_dir,i))
            shell.writelines("%s view %s/%s.sam -S -b -o %s/%s.bam\n" % \
                             (conf_dict['samtools']['path'],singel_bowtie_result_dir,i,\
                              singel_bowtie_result_dir,i))
        bam_dict[i] = os.path.join(singel_bowtie_result_dir,i,".bam")

    return bam_dict

## step3 Gene and isoform expresssion level
## based on gene expression
def gene_exp_rsem(conf_dict,samples_dict, bam_dict, result_dir):
    ## the bam file produced by bowtie2
    rsem = conf_dict['rsem']['path']
    refMrna = conf_dict['refMrna']
    rsem_shell_dir = os.path.join(result_dir, "shell/gene_exp_rsem")
    rsem_result_dir = os.path.join(result_dir, "result/rsem_result_dir")
    for i in samples_dict.keys():
        singel_rsem_shell_dir = check_dir(os.path.join(rsem_shell_dir, i))
        singel_rsem_result_dir = check_dir(os.path.join(rsem_result_dir, i))
        shell_file = os.path.join(singel_rsem_shell_dir, "work.sh")
        with open (shell_file,'w') as shell:
            shell.writelines("{bin}/rsem-calculate-expression -p 8 --bam {bam} {referMrna} {outdir}\n".format(
                bin=rsem,bam=bam_dict[i],referMrna=refMrna,outdir=singel_rsem_result_dir))
            # 统计基因表达量
            awk_cmd = "awk '{OFS=\"\\t\";if($7!=0.00)print $1,$2,$3,$5,$7}'"

            shell.writelines("{awk_cmd} {outdir}/{sample}.genes.results | grep -v '^ERCC' > {outdir}/{sample}.gene.fpkm.xls\n ".format(
                outdir=singel_rsem_result_dir,awk_cmd=awk_cmd,sample=i)
            )
            shell.writelines("{awk_cmd} {outdir}/{sample}.isoforms.results | grep -v '^ERCC' > {outdir}/{sample}.transcript.fpkm.xls\n ".format(
                outdir=singel_rsem_result_dir, awk_cmd=awk_cmd, sample=i)
            )



##step4 Deep analysis of gene expression
def gene_diff_exp(conf_dict, result_dir):
    # method 1 : PossionDis
    # PossionDis = conf_dict['PossionDis']
    # PossionDis_VS = T1&C1,T2&C2
    # PossionDis_Filter = conf_dict[PossionDis_Filter]

    # method 2 : DESeq

    # method 3 : DESeq2
    diffexp = conf_dict['diffexp']

    for key,value in diffexp['DESeq2'].items():
        print(key,value)
        assert(len(value)>=2)
    # 将上面得到的表达量的数据转化成本地矩阵文件，适合DESeq2读取
    # 执行DESeq2
    # 获取DESeq2的结果
    # 输出差异表达文件

    # method 4 : NOIseq

    # method 5 : EBsqe
def gene_pca():
    pass
def gene_co_expression_network():
    pass
def gene_cluster():
    pass



def main(args):
    if args == None:
        # parse_args()
        sys.exit(1)

    print(args)

    samples_list = args.samples
    configure_path = args.config
    output_dir = args.output


    samples = sample_parse(samples_list)

    ## 获取配置文件
    conf_dict = configure(configure_path)
    print("conf_dict:",conf_dict)

    ## 获取当前目录 为后续建立目录
    result_dir = output_dir
    ## 过滤数据
    samples_dict = filter_cutadapt(samples, conf_dict, result_dir)
    ## 基因组比对
    genomeMapping_hisat(conf_dict, samples_dict, result_dir)
    #
    # ## 基因比对
    bam_dict = gene_mapping_bowtie(conf_dict, samples_dict, result_dir)
    # ## 表达量计算
    gene_exp_rsem(conf_dict, samples_dict, bam_dict, result_dir)

    # diff
    gene_diff_exp(conf_dict)



    ##step5 GO analysis

    ##step6 KEGG analysis

if __name__ == "__main__":
    args = parse_args()
    main(args)
