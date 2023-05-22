import pandas as pd, numpy as np
import os 

configfile: "config/config.yml"


def load_blast(ff):
    df=pd.read_csv(ff,sep="\t",names=["query","subject","per_identity","alignment_length","mismatches","gap_opens","query_start","query_end","subject_start","subject_end", "evalue", "bit_score"],comment='#')
    return df

def generate_mcl_file(ff,ff_out,minbit_threshold=0.5):
    df=load_blast(ff)
    df=df.drop_duplicates(['query','subject'],keep='first') # Drop duplicate maps of same gene pair
    df['self_map']=df['query']==df['subject']

    print("Calculating minbit")
    bit_score_selfmap=df[df['self_map']][['query','bit_score']].copy()
    df=df.merge(bit_score_selfmap,left_on='query',right_on='query',how='left',suffixes=['','_selfmap_q'])
    df=df.merge(bit_score_selfmap,left_on='subject',right_on='query',how='left',suffixes=['','_selfmap_s']).drop(columns=['query_selfmap_s'])
    df['minbit']=df['bit_score']/np.minimum(df['bit_score_selfmap_q'].values,df['bit_score_selfmap_s'].values)

    mcl_input=df[df['minbit']>minbit_threshold][['query','subject','per_identity']].copy().rename(columns={'per_identity':'weight'})
    mcl_input['weight']=mcl_input['weight']/100
    mcl_input.to_csv(ff_out,sep="\t",index=False,header=False)

def load_mcl_output(ff,cluster_prefix=''):
    out=[]
    with open(ff) as f:
        for i,line in enumerate(f.read().split('\n')):
            out.extend([[cluster_prefix+str(i),_] for _ in line.split('\t')])
    return pd.DataFrame(out,columns=['cluster','gene_id'])


rule all:
    input:
        mcl_out=os.path.join(config['output'],'mcl.out')

rule all_all_blast:
    input: 
        config['input']
    output: 
        blast_out=os.path.join(config['output'],'blast.out')
    params: 
        blast_db=os.path.join(config['output'],'blast.db'),
        input_type=config['input_type']
    threads: 
        config['threads']
    conda: 'genecluster'
    shell:
        # config 
        """
        echo "Running all-all blast..."
        input_type={params.input_type}
        if [ $input_type == "faa" ]; then
            makeblastdb -in {input} -dbtype prot -out {params.blast_db}
            blastp -db {params.blast_db} -query {input} -out {output.blast_out} -num_threads {threads} -outfmt 7 -evalue 1e-5 -max_target_seqs 10000
        elif [ $input_type == "fna" ]; then
            makeblastdb -in {input} -dbtype nucl -out {params.blast_db}
            blastn -db {params.blast_db} -query {input} -out {output.blast_out} -num_threads {threads} -outfmt 7 -evalue 1e-5 -max_target_seqs 10000
        else
            echo "Input type not recognized. Please use faa or fna"
            exit 1
        fi
        """


rule prepare_mcl_file:
    input:
        blast_out=os.path.join(config['output'],'blast.out')
    output:
        mcl_in=os.path.join(config['output'],'mcl.in')
    run:
        print("Generating mcl input file...")
        generate_mcl_file(input.blast_out,output.mcl_in)


rule mcl:
    input:
        mcl_in=os.path.join(config['output'],'mcl.in')

    output:
        mcl_out=os.path.join(config['output'],'mcl.out')
    params:
        inflation=config['mcl_inflation'],
    threads:
        config['threads']
    conda: 'genecluster'
    shell:
        """
        mcl {input.mcl_in} --abc -I {params.inflation} -o {output.mcl_out} -te {threads}
        """