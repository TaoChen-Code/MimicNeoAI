# coding=<encoding name> ： # coding=utf-8
import os
import sys
import multiprocessing

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

if sys.version_info < (3, 8):
    # Python 3.8以下版本实现
    class NoDaemonPool(multiprocessing.pool.Pool):
        Process = NoDaemonProcess
else:
    # Python 3.8及以上版本实现
    class NoDaemonPool(multiprocessing.pool.Pool):
        @staticmethod
        def Process(_, *args, **kwds):
            return NoDaemonProcess(*args, **kwds)

class Pvacseq():
    def __init__(self, tool):
        self.tool = tool

    def get_hlahd_results(self,sample,output_hla):
        hlaI_set = ['A','B','C','E','F','G','H','J','K','L','V']
        hlaII_set = ['DRA','DRB1','DRB2','DRB3','DRB4','DRB5','DRB6','DRB7','DRB8','DRB9',
                     'DQA1','DQB1','DPA1','DPB1','DMA','DMB','DOA','DOB']
        
        hla_file = f"{output_hla}/{sample}/result/{sample}_final.result.txt"
        hlaI = ''
        hlaII = ''
        
        with open(hla_file, "r") as f:
            for line in f.readlines():
                new_line = line.replace("\n","").split('\t')
                
                # Process HLA-I alleles
                if new_line[0] in hlaI_set:
                    for h in new_line[1:]:
                        h_split = h.split(":")
                        if h != '-' and h != 'Not typed' and len(h_split) >= 2:
                            h_2 = ':'.join(h_split[:2])  # Keep only first two fields
                            hlaI += h_2 + ','
                
                # Process HLA-II alleles
                if new_line[0] in hlaII_set:
                    for h in new_line[1:]:
                        h_split = h.split(":")
                        if h != '-' and h != 'Not typed' and len(h_split) >= 2:
                            h_2 = ':'.join(h_split[:2])
                            # Remove HLA prefix (e.g., "DRB1" from "DRB1*15:01")
                            hlaII += h_2.split('-')[1] + ','
                            
        return hlaI, hlaII
    
    def pvacseq(self,run_sample_id, sample, output_vcf,output_hla,configure,pathes,vcf_chunk = None):
        """
        Run pVACseq pipeline on a VCF file or chunk with multiprocessing support.
        
        Parameters:
            sample (str): Sample ID/name
            output_vcf (str): Path to input VCF file
            output_hla (str): Path to HLA typing results
            configure (dict): Configuration dictionary containing:
                - args: thread count
                - others: species, mutation_calling_tool
                - path: output directory
                - step_name: step names for annotation and pvacseq
            pathes (dict): Paths to tools including:
                - path: pvacseq and iedb-install-directory
            vcf_chunk (str, optional): Specific VCF chunk to process
        
        Returns:
            str: Path to output directory
        """
        #pvacbind args
        # algoIandII = 'MHCnuggetsI MHCnuggetsII NNalign NetMHC NetMHCIIpan SMM SMMPMBEC SMMalign NetMHCpanEL NetMHCIIpanEL'
        # algoI = 'MHCnuggetsI NNalign NetMHC SMM SMMPMBEC SMMalign NetMHCpanEL'
        #algoIandII = 'NNalign NetMHC NetMHCIIpan SMM SMMPMBEC SMMalign NetMHCpanEL NetMHCIIpanEL'
        
        # # Neo algo
        # algoIandII = 'BigMHC_EL BigMHC_IM DeepImmuno MHCflurry MHCflurryEL MHCnuggetsI MHCnuggetsII NNalign \
        # NetMHC NetMHCIIpan NetMHCIIpanEL NetMHCpan NetMHCpanEL PickPocket SMM SMMPMBEC SMMalign'
        # algoI = 'NNalign NetMHC SMM SMMPMBEC SMMalign NetMHCpanEL'

        # Neo algo
        algoIandII = 'BigMHC_EL BigMHC_IM DeepImmuno MHCflurry MHCflurryEL MHCnuggetsI MHCnuggetsII NNalign \
        NetMHC NetMHCIIpan NetMHCIIpanEL NetMHCpan NetMHCpanEL PickPocket SMM SMMPMBEC'
        algoI = 'NNalign NetMHC SMM SMMPMBEC NetMHCpanEL'
        
        
        #configure
        thread = configure['args']['thread']
        species = configure['others']['species']
        output_dir = configure['path']['output_dir'] + "/" #输出文件夹
        mutation_calling_tool = configure['others']['mutation_calling_tool']
        pvacseq = pathes['path']['pvacseq']
        iedb_install_directory=pathes['path']['iedb-install-directory']
        #创建输出文件夹
        step_name_annotation = configure['step_name']['annotation']
        step_name_pvacseq = configure['step_name']['pvacseq']
        output_pvacseq = output_dir + f'/{sample}/{step_name_pvacseq}/'
        if vcf_chunk is not None:
            output_pvacseq = output_dir + f'/{sample}/{step_name_pvacseq}/chunks/chunk_{vcf_chunk}'
    
        cmd_mkdir = "mkdir -p {}".format(output_pvacseq)
        self.tool.judge_then_exec(run_sample_id,cmd_mkdir,output_pvacseq)
        #后缀
        if mutation_calling_tool == 'HaplotypeCaller':
            if species == "mouse":
                suffix = "HC.snps.indels.filtered.PASS.VEP.filtered.vcf"
            else:
                suffix = "HC.snps.indels.VQSR.PASS.VEP.filtered.vcf"
        elif mutation_calling_tool == 'Mutect2':
            suffix = "mutect.filtered.VEP.filtered.vcf"
        #read hla
        hlaI,hlaII = self.get_hlahd_results(sample,output_hla) 
        #HLA = 'HLA-A*31:01,HLA-A*31:01,HLA-B*39:01,HLA-C*07:02,DQB1*06:01,DRB1*08:03'
        if hlaI != '' and hlaII != '':
            HLA = hlaI + hlaII[:-1]
    
            # Build the pvacseq command
            input_vcf = f"{output_dir}/{sample}/{step_name_annotation}/{sample}.{suffix}"
            running_thread = thread
            if vcf_chunk is not None:
                input_vcf = f"{output_dir}/{sample}/{step_name_annotation}/chunks/{sample}.{suffix[:-4]}.chunk_{vcf_chunk}.vcf"
                running_thread = 1
            
            cmd = f"{pvacseq} run {input_vcf} {sample} {HLA} {algoIandII} {output_pvacseq} -e1 8,9,10 -e2 15 --iedb-install-directory {iedb_install_directory} -t {running_thread} --fasta-size 100000"
            # cmd_export1 = "export TF_CPP_MIN_LOG_LEVEL=2"
            # cmd_export2 = "export CUDA_VISIBLE_DEVICES=\"\""
            # self.tool.exec_cmd(cmd_export1,sample)
            # self.tool.exec_cmd(cmd_export2,sample)
            self.tool.judge_then_exec(run_sample_id,cmd,f"{output_pvacseq}/combined/{sample}.all_epitopes.tsv",vcf_chunk)
        else:
            self.tool.write_log("hlahd results is all Not_typed!","error")
        return output_pvacseq
          
    def run_pvacseq_parallel(self,run_sample_id,sample,output_vcf,output_hla,configure,pathes):
        """
        Execute pvacseq in parallel across 10 VCF chunks.
        
        Parameters:
            sample (str): Sample ID/name
            output_vcf (str): Path to input VCF file
            output_hla (str): Path to HLA typing results
            configure (dict): Configuration dictionary
            pathes (dict): Paths to tools
        """
        thread = configure['args']['thread']
        output_dir = configure['path']['output_dir'] + "/" #输出文件夹
        step_name_annotation = configure['step_name']['annotation']
        step_name_pvacseq = configure['step_name']['pvacseq']
        output_pvacseq = output_dir + f'/{sample}/{step_name_pvacseq}/'
        pool = NoDaemonPool(thread)
        for i in range(thread):
            chunk = str(i)
            pool.apply_async(self.pvacseq,(run_sample_id, sample, output_vcf, output_hla, configure, pathes, chunk,),error_callback=self.tool.print_pool_error)
        pool.close()
        pool.join()
        
        # merge
        self.merge_chunk_results(sample, output_pvacseq, thread)

    def merge_chunk_results(self, sample, output_dir, num_chunks):
        """
        Merge {sample}.all_epitopes.tsv from available chunks into a single file.

        - Only merges files that exist and contain data rows (beyond the header).
        - Writes the header exactly once.
        - Logs a summary; if no usable files are found, logs an error.

        Args:
            sample (str): Sample ID.
            output_dir (str): Base output directory.
            num_chunks (int): Number of chunks to check.
        """
        merged_dir = f"{output_dir}/combined"
        os.makedirs(merged_dir, exist_ok=True)

        merged_epitopes = f"{merged_dir}/{sample}.merged.all_epitopes.tsv"

        # Candidate paths for each chunk
        candidates = [
            f"{output_dir}/chunks/chunk_{i}/combined/{sample}.all_epitopes.tsv"
            for i in range(num_chunks)
        ]

        # Keep only existing, non-empty files
        existing = [p for p in candidates if os.path.exists(p) and os.path.getsize(p) > 0]

        if not existing:
            self.tool.write_log("No chunk result files found to merge.", "error")
            return

        wrote_header = False
        total_rows = 0
        used_files = []

        try:
            with open(merged_epitopes, "w") as out_f:
                for path in existing:
                    with open(path, "r") as in_f:
                        header = in_f.readline()
                        body_lines = in_f.readlines()

                        # Skip files that only contain a header
                        if not body_lines:
                            continue

                        if not wrote_header:
                            out_f.write(header)
                            wrote_header = True

                        out_f.writelines(body_lines)
                        total_rows += len(body_lines)
                        used_files.append(path)

            if total_rows == 0:
                self.tool.write_log("All chunk files contain only headers (no data).", "error")
                try:
                    os.remove(merged_epitopes)
                except OSError:
                    pass
                return

            self.tool.write_log(
                f"Merged {len(used_files)} files with {total_rows} data rows -> {merged_epitopes}",
                "info"
            )
        except IOError as e:
            raise Exception(f"Error while merging files: {e}")


  
