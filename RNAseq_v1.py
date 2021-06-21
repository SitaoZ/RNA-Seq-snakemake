#coding:utf-8
#!/usr/bin/python

import os,sys
from optparse import OptionParser

##python 2.7版本支持
# for python 3.6 
parser = OptionParser(usage="%prog [options] Samples_List Configure")

parser.add_option("-t", "--timeout",
                action = "store",
                type = 'int',
                dest = "timeout",
                default = None,
                help="Specify annalysis execution time limit"
                )
parser.add_option("-u", "--url",
                action = "store_true",
                dest = "url",
                default = False,
                help = "Specify if the target is an URL"
                )

(options,args) = parser.parse_args()

def Sample(Sample_Path):
    Samples = {}
    with open(Sample_Path, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if not len(line) or line.startswith('#'):
                continue
            else:
                SS = line.split('\t')
                Samples[SS[0]] = [SS[1], SS[2]]
    return Samples

def Configure(Configure_Path):
    import re
    Conf_Dict = {}
    Conf = open(Configure_Path,'r')
    for line in Conf.readlines():
        line = line.strip()
        if not len(line) or line.startswith('#'):
            continue
        else:
            pattern = re.compile(r'\S+')
            Array = line.split('=')
            match1 = pattern.search(Array[0])
            match2 = pattern.search(Array[1])
            if match1 and match2 :
                Conf_Dict[match1.group(0)] = match2.group(0)
    return Conf_Dict

def Check_dir(path):
    import os
    if not os.path.exists(path):
        os.makedirs(path)
        return path
    else:
        print path,"exists or file_name unusual!"


# def Make_Dir(pwd):
#
#     ## layer1 directory
#     Check_dir(pwd+"/"+"result")
#     Check_dir(pwd+"/"+"shell")
#
#     ## layer 2 directory
#     Check_dir(pwd + "/result/Filter")
#     Check_dir(pwd + "/result/GenomeMapping_HISAT")
#     Check_dir(pwd + "/result/GeneMapping_Bowtie")
#     Check_dir(pwd + "/result/GeneExp_RSEM")
#
#
#     Check_dir(pwd + "/shell/Filter")
#     Check_dir(pwd + "/shell/GenomeMapping_HISAT")
#     Check_dir(pwd + "/shell/GeneMapping_Bowtie")
#     Check_dir(pwd + "/shell/GeneExp_RSEM")




##step1 filter
def Filter_SOAPnuke(Samples,Conf_Dict):
    Samples_Dict = {}
    Filter_Shell_Dir = pwd + "/shell/Filter_SOAPnuke"
    Filter_Result_Dir = pwd + "/result/Filter_SOAPnuke"

    for i in Samples.keys():
        Singel_Filter_Shell_dir = Check_dir(Filter_Shell_Dir+"/"+i)
        Singel_Filter_Result_dir = Check_dir(Filter_Result_Dir+"/"+i)

        Shell_File = Singel_Filter_Shell_dir+"/work.sh"
        with open(Shell_File,'w') as shell:
            shell.writelines("zcat %s %s |gzip > %s/%s.fq.gz\n"%(Samples[i][0],Samples[i][1],Singel_Filter_Result_dir,i))
            shell.writelines("%s %s %s/%s.fq.gz\n" % (Conf_Dict['SOAPnuke'],Conf_Dict['SOAPnuke_Options'],Singel_Filter_Result_dir,i))
        Samples_Dict[i] = "%s/%s.fq.gz" % (Singel_Filter_Result_dir,i)
    return Samples_Dict


##step2 alignment
# hisat alignment of genome reference
def GenomeMapping_HISAT(Conf_Dict,Samples_Dict):
    HISAT = Conf_Dict['HISAT']
    HISAT_options_PE = Conf_Dict['HISAT_options_PE']
    HISAT_Index = Conf_Dict['HISAT_Index']

    HISAT_Shell_Dir = pwd + "/shell/HISAT"
    HISAT_Result_Dir = pwd + "/result/HISAT"

    for i in Samples_Dict.keys():
        Singel_HISAT_Shell_dir = Check_dir(HISAT_Shell_Dir + "/" + i)
        Singel_HISAT_Result_dir = Check_dir(HISAT_Result_Dir + "/" + i)

        Shell_File = Singel_HISAT_Shell_dir + "/work.sh"
        with open(Shell_File, 'w') as shell:
            shell.writelines("%s %s -x %s -U %s -S %s.sam 2>%s.stat\n"%(HISAT,HISAT_options_PE,HISAT_Index,Samples_Dict[i],Singel_HISAT_Result_dir,Singel_HISAT_Result_dir))
            shell.writelines("%s view %s.sam -b -S -o %s.bam\n"%(Conf_Dict['samtools'],Singel_HISAT_Result_dir,Singel_HISAT_Result_dir))



# Bowtie alignment of gene reference (for gene expression analysis)
def GeneMapping_Bowtie(Conf_dict,Samples_Dict):
    Bam_Dict = {}
    Bowtie = Conf_Dict['Bowtie']
    GeneBowtie2Index = Conf_Dict['Bowtie2Index']
    GeneBowtie2_Options_PE = Conf_Dict['GeneBowtie2_Options_PE']

    Bowtie_Shell_Dir = pwd + "/shell/GeneMapping_Bowtie"
    Bowtie_Result_Dir = pwd + "/result/GeneMapping_Bowtie"
    for i in Samples_Dict.keys():
        Singel_Bowtie_Shell_Dir = Check_dir(Bowtie_Shell_Dir + "/" + i)
        Singel_Bowtie_Result_Dir = Check_dir(Bowtie_Result_Dir + "/" + i)
        Shell_File = Singel_Bowtie_Shell_Dir + "/work.sh"
        with open(Shell_File,'w') as shell:
            shell.writelines("%s %s -x %s -U %s > %s/%s.sam\n"%(Bowtie,GeneBowtie2_Options_PE,Conf_Dict['Bowtie2Index'],Samples_Dict[i],Singel_Bowtie_Result_Dir,i))
            shell.writelines("%s view %s/%s.sam -S -b -o %s/%s.bam\n"%(Conf_Dict['samtools'],Singel_Bowtie_Result_Dir,i,Singel_Bowtie_Result_Dir,i))
        Bam_Dict[i] = Singel_Bowtie_Result_Dir+"/"+i+".bam"
    return Bam_Dict

##step3 Gene and isoform expresssion level
## based on gene expression
def GeneExp_RSEM(Conf_Dict,Samples_Dict,Bam_Dict):
    ## the bam file produced by bowtie2
    RSEM = Conf_Dict['RSEM']
    RSEM_Shell_Dir = pwd + "/shell/GeneExp_RSEM"
    RSEM_Result_Dir = pwd + "/result/RSEM_Result_Dir"
    for i in Samples_Dict.keys():
        Singel_RSEM_Shell_Dir = Check_dir(RSEM_Shell_Dir+"/"+i)
        Singel_RSEM_Result_Dir = Check_dir(RSEM_Result_Dir + "/" + i)
        Shell_File = Singel_RSEM_Shell_Dir+"/work.sh"
        with open (Shell_File,'w') as shell:
            shell.writelines("%s --bam %s -p 8 %s %s\n"%(RSEM,Bam_Dict[i],Conf_Dict['Gene'],Singel_RSEM_Result_Dir))
            # 统计基因表达量
            shell.writelines("")




##step4 Deep analysis of gene expression
def DeneDiffExp(Conf_Dict):

    # method 1 : PossionDis
    PossionDis = Conf_Dict['PossionDis']
    PossionDis_VS = T1&C1,T2&C2
    PossionDis_Filter = Conf_Dict[PossionDis_Filter]

    # method 2 : DESeq

    # method 3 : DESeq2

    # method 4 : NOIseq

    # method 5 : EBsqe
def Gene_PCA():
    pass
def Gene_Co_Expression_Network():
    pass
def Gene_Cluster():
    pass






##step5 GO analysis

##step6 KEGG analysis

if __name__ == "__main__":
    #(options, args) = parser.parse_args()
    print '第一个位置参数',args[0]
    print '第二个位置参数',args[1]
    Samples_List = args[0]
    Configure_Path = args[1]

    Samples = Sample(Samples_List)

## 获取配置文件
    Conf_Dict = Configure(Configure_Path)

## 获取当前目录 为后续建立目录
    pwd = os.getcwd()

## 过滤数据
    Samples_Dict = Filter_SOAPnuke(Samples,Conf_Dict)
## 基因组比对
    GenomeMapping_HISAT(Conf_Dict,Samples_Dict)

## 基因比对
    Bam_Dict = GeneMapping_Bowtie(Conf_Dict,Samples_Dict)
## 表达量计算
    GeneExp_RSEM(Conf_Dict,Samples_Dict,Bam_Dict)
