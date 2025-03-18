def get_data_for_pvacseq(blast_file,pvacbind_file,output_blastx):
    context = ""
    with open(f"{output_blastx}/{blast_file}",'r') as f:
        lines = f.readlines()
        for i,line in enumerate(lines):
            id = f">index:{i}|" + line.split("\t")[0]
            peptide = line.split("\t")[2]
            if peptide != '':
                context += id + '\n' + peptide + '\n'
            #break
    context = context[:-1]

    with open(f"{output_blastx}/{pvacbind_file}","w") as f:
        f.write(context)