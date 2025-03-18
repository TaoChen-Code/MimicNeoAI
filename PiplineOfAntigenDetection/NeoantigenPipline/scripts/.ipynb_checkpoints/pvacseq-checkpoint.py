# coding=<encoding name> ： # coding=utf-8
#mouse
import os
import pandas as pd

def get_hlahd_results(sample,output_hla):
  hlaI_set = ['A','B','C','E','F','G','H','J','K','L','V']
  hlaII_set = ['DRA','DRB1','DRB2','DRB3','DRB4','DRB5','DRB6','DRB7','DRB8','DRB9','DQA1','DQB1','DPA1','DPB1','DMA','DMB','DOA','DOB']
  hla_file = f"{output_hla}/{sample}/result/{sample}_final.result.txt"
  hlaI = ''
  hlaII = ''
  with open(hla_file,"r") as f:
    for line in f.readlines():
      new_line = line.replace("\n","").split('\t')
      #hlaI
      if new_line[0] in hlaI_set:
        for h in new_line[1:]:
          h_split = h.split(":")
          if h != '-' and h != 'Not typed' and len(h_split) >= 2:
            h_2 = ':'.join(h_split[:2])  # h_2为只保留两位的hla位点，因为pvacbind不支持三位
            hlaI += h_2 + ','
      if new_line[0] in hlaII_set:
        for h in new_line[1:]:
          h_split = h.split(":")
          if h != '-' and h != 'Not typed' and len(h_split) >= 2:
            h_2 = ':'.join(h_split[:2])
            hlaII += h_2.split('-')[1] + ','
  return hlaI,hlaII

def pvacseq(sample,output_vcf,output_hla,configure,pathes,tool):
    #pvacbind args
    # algoIandII = 'MHCnuggetsI MHCnuggetsII NNalign NetMHC NetMHCIIpan SMM SMMPMBEC SMMalign NetMHCpanEL NetMHCIIpanEL'
    #algoI = 'MHCnuggetsI NNalign NetMHC SMM SMMPMBEC SMMalign NetMHCpanEL'
    algoIandII = 'NNalign NetMHC NetMHCIIpan SMM SMMPMBEC SMMalign NetMHCpanEL NetMHCIIpanEL'
    algoI = 'NNalign NetMHC SMM SMMPMBEC SMMalign NetMHCpanEL'
    #configure
    thread = configure['args']['thread']
    species = configure['others']['species']
    sif_file = pathes['path']['sif_file']
    output_dir = configure['path']['output_dir'] + "/" #输出文件夹
    scp_dir = configure['path']['scp_dir']
    mutation_calling_tool = configure['others']['mutation_calling_tool']
    pvacseq = pathes['path']['pvacseq']
    iedb_install_directory=pathes['path']['iedb-install-directory']
    #创建输出文件夹
    output_pvacseq = output_dir + '/{}/10-pvacseq/'.format(sample)
    # cmd_mkdir = "mkdir {}".format(output_pvacseq)
    #tool.judge_then_exec(sample,cmd_mkdir,output_pvacseq)
    #后缀
    if mutation_calling_tool == 'HaplotypeCaller':
        if species == "mouse":
            suffix = "HC.snps.indels.filtered.PASS.VEP.filtered.vcf"
        else:
            suffix = "HC.snps.indels.VQSR.PASS.VEP.filtered.vcf"
    elif mutation_calling_tool == 'Mutect2':
        suffix = "mutect.filtered.VEP.filtered.vcf"
    #read hla
    hlaI,hlaII = get_hlahd_results(sample,output_hla) 
    #hlaII
    #HLA = 'HLA-A*31:01,HLA-A*31:01,HLA-B*39:01,HLA-C*07:02,DQB1*06:01,DRB1*08:03'
    if hlaI != '' and hlaII != '':
        HLA = hlaI + hlaII[:-1]
        #pvacseq cmd
        cmd = f"{pvacseq} run {output_dir}/{sample}/08-vep/{sample}.{suffix} {sample} {HLA} {algoIandII} {output_pvacseq} -e1 8,9,10 -e2 15 --iedb-install-directory {iedb_install_directory} -t {thread}"
        pvacseq_sh = tool.write2shell(cmd,"pvacseq")
        tool.judge_then_exec(sample,cmd,f"{output_pvacseq}/combined/{sample}.filtered.tsv")
        
        #scp
        # tool.scp2ifdgm(f"{output_vcf}/{sample}.{suffix}",f"{scp_dir}/08.vep/",sample)
        # tool.scp2ifdgm(pvacseq_sh,scp_dir,sample)

        #docker cmd
        # cmd_docker = f"docker run -P -v {output_vcf}:/data_dir/ -v {output_pvacseq}:/output_file/ griffithlab/pvactools pvacseq run /data_dir/{sample}.{suffix} {sample} {HLA} {algoIandII} /output_file/ -e1 9 -e2 15 --iedb-install-directory /opt/iedb -t 5"
        # tool.write2shell(cmd_docker,"docker")
    else:
        tool.write_log("hlahd results is all Not_typed!","error")
    return output_pvacseq
      
def pvacbind_for_mouse(sample,fa_dir,output_dir,configure,pathes,tool):
  pass

  
  


  
