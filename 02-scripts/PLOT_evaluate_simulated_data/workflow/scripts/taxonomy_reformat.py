
#########################################################
# 1. Importing libraries
#########################################################

import polars as pl
import polars.selectors as cs
import pandas as pd
from pathlib import Path
import pytaxonkit
import sys

#########################################################
# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")


#########################################################
# 2. Functions
#########################################################

def parse_output_metabuli_files(metabuli_file, big_unwanted_list):
    """
    Parse output centrifuge files and extract taxonomy and species ID from read.
    
    Parameters:
    metabuli_files (list): List of FASTQ file paths.
    
    Returns:
    pd.DataFrame: DataFrame with columns for file name, read ID, taxonomy ID, and species ID.
    """
    
    
    if not metabuli_file.with_suffix('.output.parquet').exists():
        metabuli_df = pl.read_csv(metabuli_file,
                            separator="\t",
                            has_header = False,
                            new_columns = ['Classified', 'Read ID', 'Taxonomy ID', 'Read Length', 'ID score', 'Rank', 'List taxid:kmer'],
                            dtypes = [pl.Int32, pl.Utf8, pl.Int32, pl.Int32, pl.Float32, pl.Utf8, pl.Utf8],
                            ).with_columns(
                                [
                                    pl.col("Read ID").str.split(':').list.first().alias("Species ID"),
                                    pl.lit(metabuli_file.stem.split('.')[0]).alias("File Name"),
                                    pl.col("Classified").replace_strict([0,1], ['U', 'C']),
                                ]
                            )

        metabuli_df = metabuli_df.filter(
            ~pl.col('Species ID').is_in(big_unwanted_list)
        )

        # metabuli_df.to_pandas().to_csv(metabuli_file.with_suffix('.output.tsv.gz'), index=False, sep='\t')
        metabuli_df.write_parquet(metabuli_file.with_suffix('.output.parquet'))
    else:
        metabuli_df = pl.read_parquet(metabuli_file.with_suffix('.output.parquet'))

    return metabuli_df

#########################################################

def parse_output_centrifuge_files(centrifuge_file, big_unwanted_list):
    """
    Parse output centrifuge files and extract taxonomy and species ID from read.
    
    Parameters:
    centrifuge_files (list): List of FASTQ file paths.
    
    Returns:
    pd.DataFrame: DataFrame with columns for file name, read ID, taxonomy ID, and species ID.
    """
    
    
    if not centrifuge_file.with_suffix('.parquet').exists():
        centrifuge_df = pl.read_csv(centrifuge_file,
                            separator="\t",
                            has_header = True,
                            dtypes = [pl.Utf8, pl.Utf8, pl.Int32, pl.Int32, pl.Int32, pl.Int32, pl.Int32, pl.Int32],
                            ).rename(
                                {
                                    "readID": "Read ID",
                                    "taxID": "Taxonomy ID",
                                    "seqID": "Sequence ID",   
                                }
                            ).with_columns(
                                [
                                    pl.col("Read ID").str.split(':').list.first().alias("Species ID"),
                                    pl.lit(centrifuge_file.stem.split('.')[0]).alias("File Name"),
                                ]
                            )

        centrifuge_df = centrifuge_df.filter(
            ~pl.col('Species ID').is_in(big_unwanted_list)
        )

        # centrifuge_df.to_pandas().to_csv(centrifuge_file.with_suffix('.tsv.gz'), index=False, sep='\t')
        centrifuge_df.write_parquet(centrifuge_file.with_suffix('.parquet'))
    else:
        centrifuge_df = pl.read_parquet(centrifuge_file.with_suffix('.parquet'))

    return centrifuge_df
    

#########################################################

def parse_output_kraken2_files(kraken2_file, big_unwanted_list):
    """
    Parse output kraken2 files and extract taxonomy and species ID from read.
    
    Parameters:
    kraken2_files (list): List of FASTQ file paths.
    
    Returns:
    pd.DataFrame: DataFrame with columns for file name, read ID, taxonomy ID, and species ID.
    """
    
    
    if not kraken2_file.with_suffix('.parquet').exists():
        kraken_df = pl.read_csv(kraken2_file, 
                            separator="\t",
                            has_header = False,
                            new_columns = ['Classified', 'Read ID', 'Taxonomy ID', 'Read length', 'LCA kmer'],
                            dtypes = [pl.Utf8, pl.Utf8, pl.Utf8, pl.Utf8, pl.Utf8],
                            ).with_columns(
                                [
                                    pl.col("Read ID").str.split(':').list.first().alias("Species ID"),
                                    pl.col("Taxonomy ID").str.extract(r'taxid ([0-9]+)', 1),
                                    pl.lit(kraken2_file.stem.split('.')[0]).alias("File Name"),
                                ]
                            )

        kraken_df = kraken_df.filter(
            ~pl.col('Species ID').is_in(big_unwanted_list)
        )

        # kraken_df.to_pandas().to_csv(kraken2_file.with_suffix('.tsv.gz'), index=False, sep='\t')
        kraken_df.write_parquet(kraken2_file.with_suffix('.parquet'))
    else:
        kraken_df = pl.read_parquet(kraken2_file.with_suffix('.parquet'))

    return kraken_df

#########################################################

def parse_output_krakenuniq_files(krakenuniq_file, big_unwanted_list):
    """
    Parse output kraken2 files and extract taxonomy and species ID from read.
    
    Parameters:
    krakenuniq_files (list): List of FASTQ file paths.
    
    Returns:
    pd.DataFrame: DataFrame with columns for file name, read ID, taxonomy ID, and species ID.
    """
    
    
    if not krakenuniq_file.with_suffix('.parquet').exists():
        kraken_df = pl.read_csv(krakenuniq_file, 
                            separator="\t",
                            has_header = False,
                            new_columns = ['Classified', 'Read ID', 'Taxonomy ID', 'Read length', 'LCA kmer'],
                            dtypes = [pl.Utf8, pl.Utf8, pl.Int32, pl.Int32, pl.Utf8],
                            ).with_columns(
                                [
                                    pl.col("Read ID").str.split(':').list.first().alias("Species ID"),
                                    pl.lit(krakenuniq_file.stem.split('.')[0]).alias("File Name"),
                                ]
                            )

        kraken_df = kraken_df.filter(
            ~pl.col('Species ID').is_in(big_unwanted_list)
        )

        if snakemake.params.seq2name != '':
            seq2name = pd.read_table(
                snakemake.params.seq2name,
                header=None,
                names=['tmp_ID', 'Taxonomy ID']
            ).set_index('Taxonomy ID')['tmp_ID'].to_dict()

            seq2name_old = pd.read_table(
                snakemake.params.seq2name_old,
                header=None,
                names=['tmp_ID', 'Taxonomy ID']
            ).set_index('tmp_ID')['Taxonomy ID'].to_dict()

            kraken_df = kraken_df.with_columns(
                [
                    pl.col("Taxonomy ID").replace_strict(seq2name, 
                        default=pl.col("Taxonomy ID")).alias("Tmp Species ID"),
                ]
            )

            kraken_df = kraken_df.with_columns(
                [
                    pl.col("Tmp Species ID").replace_strict(seq2name_old,
                        default=pl.col("Taxonomy ID")).alias("Tmp Taxonomy ID"),
                ]
            )

            kraken_df = kraken_df.drop(
                [
                    "Taxonomy ID",
                    "Tmp Species ID"
                ]
            )

            kraken_df = kraken_df.rename(
                {
                    "Tmp Taxonomy ID": "Taxonomy ID"
                }
            )

        # kraken_df.to_pandas().to_csv(krakenuniq_file.with_suffix('.tsv.gz'), index=False, sep='\t')
        kraken_df.write_parquet(krakenuniq_file.with_suffix('.parquet'))
    else:
        kraken_df = pl.read_parquet(krakenuniq_file.with_suffix('.parquet'))

    return kraken_df

#########################################################

def transform_taxo_centrifuge(df):
    """
    Transform the taxo DataFrame
    
    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to be transformed
    
    Returns
    -------
    pl.DataFrame
        The transformed DataFrame
    """
    
    zero_ranks = 'superkingdom;clade;kingdom;phylum;class;order;family;subfamily;genus;species;no rank'
    zero_lineage = 'Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined'
    # sopEPhi = 'Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Peduoviridae;Felsduovirus;Felsduovirus SopEphi'
    # sopEPhi_rank = 'superkingdom;clade;kingdom;phylum;class;family;genus;species'
    # medusavirus = 'Viruses;Varidnaviria;Bamfordvirae;Nucleocytoviricota;Megaviricetes;Mamonoviridae;Medusavirus;Medusavirus medusae'
    # medusavirus_rank = 'superkingdom;clade;kingdom;phylum;class;family;genus;species'
    # mischivirus = 'Viruses;Riboviria;Orthornavirae;Pisuviricota;Pisoniviricetes;Picornavirales;Picornaviridae;Caphthovirinae;Mischivirus;Mischivirus E'
    # mischivirus_rank = 'superkingdom;clade;kingdom;phylum;class;order;family;subfamily;genus;species'
    # nudivirus = 'Viruses;Naldaviricetes;Lefavirales;Nudiviridae;Betanudivirus;Betanudivirus hezeae'
    # nudivirus_rank = 'superkingdom;clade;order;family;genus;species'
    # sicinivirus = 'Viruses;Riboviria;Orthornavirae;Pisuviricota;Pisoniviricetes;Picornavirales;Picornaviridae;Kodimesavirinae;Sicinivirus;Sicinivirus A'
    # sicinivirus_rank = 'superkingdom;clade;kingdom;phylum;class;order;family;subfamily;genus;species'
    # Maribacter = 'Bacteria;Bacteroidota;Flavobacteriia;Flavobacteriales;Flavobacteriaceae;Maribacter;unclassified Maribacter'
    # Maribacter_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Dyella = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Lysobacterales;Rhodanobacteraceae;Dyella;unclassified Dyella'
    # Dyella_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Limnochorda = 'Bacteria;Bacillota;Limnochordia;Limnochordales;Limnochordaceae;Limnochorda;unclassified Limnochorda'
    # Limnochorda_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Lentzea = 'Bacteria;Actinomycetota;Actinobacteria;Actinomycetes;Pseudonocardiales;Pseudonocardiaceae;Lentzea;unclassified Lentzea'
    # Lentzea_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Vagococcus = 'Bacteria;Bacillota;Bacilli;Lactobacillales;Enterococcaceae;Vagococcus;unclassified Vagococcus'
    # Vagococcus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Chryseobacterium = 'Bacteria;Bacteroidota;Flavobacteriia;Flavobacteriales;Weeksellaceae;Chryseobacterium;unclassified Chryseobacterium'
    # Chryseobacterium_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Sedimentibacter = 'Bacteria;Bacillota;Tissierellia;Sedimentibacteriales;Sedimentibacter;unclassified Sedimentibacter'
    # Sedimentibacter_rank = 'superkingdom;phylum;class;genus;species'
    # Paenibacillus_silvisoli = 'Bacteria;Bacillota;Bacilli;Bacillales;Paenibacillaceae;Paenibacillus;Paenibacillus silvisoli'
    # Paenibacillus_silvisoli_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Haloarcula = 'Archaea;Euryarchaeota;Halobacteria;Halobacteriales;Halobacteriaceae;Haloarcula;unclassified Haloarcula'
    # Haloarcula_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Lelliottia = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Lelliottia;unclassified Lelliottia'
    # Lelliottia_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Pseudomonas = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;unclassified Pseudomonas'
    # Pseudomonas_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Klebsiella = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Klebsiella;unclassified Klebsiella'
    # Klebsiella_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Pseudoalteromonas = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Alteromonadales;Pseudoalteromonadaceae;Pseudoalteromonas;unclassified Pseudoalteromonas'
    # Pseudoalteromonas_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Sinanaerobacter = 'Bacteria;Bacillota;Clostridia;Eubacteriales;Sinanaerobacter;unclassified Sinanaerobacter'
    # Sinanaerobacter_rank = 'superkingdom;phylum;class;order;genus;species'
    # Streptomyces = 'Bacteria;Actinomycetota;Actinomycetes;Kitasatosporales;Streptomycetaceae;Streptomyces;unclassified Streptomyces'
    # Streptomyces_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Thermococcus = 'Archaea;Euryarchaeota;Thermococci;Thermococcales;Thermococcaceae;Thermococcus;unclassified Thermococcus'
    # Thermococcus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Bacteroides = 'Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;unclassified Bacteroides'
    # Bacteroides_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Qipengyuania = 'Bacteria;Pseudomonadota;Alphaproteobacteria;Sphingomonadales;Erythrobacteraceae;Qipengyuania;unclassified Qipengyuania'
    # Qipengyuania_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Nguyenibacter = 'Bacteria;Pseudomonadota;Alphaproteobacteria;Rhodospirillales;Acetobacteraceae;Nguyenibacter;unclassified Nguyenibacter'
    # Nguyenibacter_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Grimontia = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Vibrionales;Vibrionaceae;Grimontia;unclassified Grimontia'
    # Grimontia_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Novosphingobium = 'Bacteria;Pseudomonadota;Alphaproteobacteria;Sphingomonadales;Sphingomonadaceae;Novosphingobium;unclassified Novosphingobium'
    # Novosphingobium_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Akkermansia = 'Bacteria;Verrucomicrobiota;Verrucomicrobiia;Verrucomicrobiales;Akkermansiaceae;Akkermansia;unclassified Akkermansia'
    # Akkermansia_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Pseudarthrobacter = 'Bacteria;Actinomycetota;Actinomycetes;Micrococcales;Micrococcaceae;Pseudarthrobacter;unclassified Pseudarthrobacter'
    # Pseudarthrobacter_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Halalkalicoccus = 'Archaea;Euryarchaeota;Halobacteria;Halobacteriales;Halococcaceae;Halalkalicoccus;unclassified Halalkalicoccus'
    # Halalkalicoccus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Mycoplasma = 'Bacteria;Mycoplasmota;Mollicutes;Mycoplasmatales;Mycoplasmataceae;Mycoplasma;unclassified Mycoplasma'
    # Mycoplasma_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Rhodococcus = 'Bacteria;Actinomycetota;Actinomycetes;Myobacteriales;Nocardiaceae;Rhodococcus;unclassified Rhodococcus'
    # Rhodococcus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Thiobacillus = 'Bacteria;Pseudomonadota;Betaproteobacteria;Nitrosomonadales;Thiobacillaceae;Thiobacillus;unclassified Thiobacillus'
    # Thiobacillus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Caproicibacter = 'Bacteria;Bacillota;Clostridia;Eubacteriales;Acutalibacteraceae;Caproicibacter;unclassified Caproicibacter'
    # Caproicibacter_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Xanthomonas = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Lysobacterales;Lysobacteraceae;Xanthomonas;unclassified Xanthomonas'
    # Xanthomonas_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Bacillus = 'Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;unclassified Bacillus'
    # Bacillus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Serratia = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Yersiniaceae;Serratia;unclassified Serratia'
    # Serratia_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Exiguobacterium = 'Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Exiguobacterium;unclassified Exiguobacterium'
    # Exiguobacterium_rank = 'superkingdom;phylum;class;order;family;genus;species'

    df.loc[df.TaxID == 0, ['FullLineage', 'FullLineageRanks']] = [zero_lineage, zero_ranks]
    df.loc[df.TaxID == 1, ['FullLineage', 'FullLineageRanks']] = [zero_lineage, zero_ranks]
    # df.loc[df.TaxID == 3118720, ['FullLineage', 'FullLineageRanks']] = [sopEPhi, sopEPhi_rank]
    # df.loc[df.TaxID == 3114988, ['FullLineage', 'FullLineageRanks']] = [medusavirus, medusavirus_rank]
    # df.loc[df.TaxID == 2870396, ['FullLineage', 'FullLineageRanks']] = [mischivirus, mischivirus_rank]
    # df.loc[df.TaxID == 3116536, ['FullLineage', 'FullLineageRanks']] = [nudivirus, nudivirus_rank]
    # df.loc[df.TaxID == 2848078, ['FullLineage', 'FullLineageRanks']] = [sicinivirus, sicinivirus_rank]
    # df.loc[df.TaxID == 3053613, ['FullLineage', 'FullLineageRanks']] = [Maribacter, Maribacter_rank]
    # df.loc[df.TaxID == 3069105, ['FullLineage', 'FullLineageRanks']] = [Dyella, Dyella_rank]
    # df.loc[df.TaxID == 3109564, ['FullLineage', 'FullLineageRanks']] = [Limnochorda, Limnochorda_rank]
    # df.loc[df.TaxID == 3108822, ['FullLineage', 'FullLineageRanks']] = [Lentzea, Lentzea_rank]
    # df.loc[df.TaxID == 3109030, ['FullLineage', 'FullLineageRanks']] = [Vagococcus, Vagococcus_rank]
    # df.loc[df.TaxID == 3116594, ['FullLineage', 'FullLineageRanks']] = [Chryseobacterium, Chryseobacterium_rank]
    # df.loc[df.TaxID == 3109366, ['FullLineage', 'FullLineageRanks']] = [Sedimentibacter, Sedimentibacter_rank]
    # df.loc[df.TaxID == 3110539, ['FullLineage', 'FullLineageRanks']] = [Paenibacillus_silvisoli, Paenibacillus_silvisoli_rank]
    # df.loc[df.TaxID == 3109565, ['FullLineage', 'FullLineageRanks']] = [Limnochorda, Limnochorda_rank]
    # df.loc[df.TaxID == 3111776, ['FullLineage', 'FullLineageRanks']] = [Haloarcula, Haloarcula_rank]
    # df.loc[df.TaxID == 3110110, ['FullLineage', 'FullLineageRanks']] = [Lelliottia, Lelliottia_rank]
    # df.loc[df.TaxID == 3110772, ['FullLineage', 'FullLineageRanks']] = [Pseudomonas, Pseudomonas_rank]
    # df.loc[df.TaxID == 3111629, ['FullLineage', 'FullLineageRanks']] = [Klebsiella, Klebsiella_rank]
    # df.loc[df.TaxID == 3110111, ['FullLineage', 'FullLineageRanks']] = [Pseudomonas, Pseudomonas_rank]
    # df.loc[df.TaxID == 3112573, ['FullLineage', 'FullLineageRanks']] = [Pseudoalteromonas, Pseudoalteromonas_rank]
    # df.loc[df.TaxID == 3111540, ['FullLineage', 'FullLineageRanks']] = [Sinanaerobacter, Sinanaerobacter_rank]
    # df.loc[df.TaxID == 3111774, ['FullLineage', 'FullLineageRanks']] = [Streptomyces, Streptomyces_rank]
    # df.loc[df.TaxID == 3111325, ['FullLineage', 'FullLineageRanks']] = [Thermococcus, Thermococcus_rank]
    # df.loc[df.TaxID == 3117552, ['FullLineage', 'FullLineageRanks']] = [Bacteroides, Bacteroides_rank]
    # df.loc[df.TaxID == 3113984, ['FullLineage', 'FullLineageRanks']] = [Qipengyuania, Qipengyuania_rank]
    # df.loc[df.TaxID == 3049350, ['FullLineage', 'FullLineageRanks']] = [Nguyenibacter, Nguyenibacter_rank]
    # df.loc[df.TaxID == 3111011, ['FullLineage', 'FullLineageRanks']] = [Grimontia, Grimontia_rank]
    # df.loc[df.TaxID == 3109595, ['FullLineage', 'FullLineageRanks']] = [Novosphingobium, Novosphingobium_rank]
    # df.loc[df.TaxID == 3115152, ['FullLineage', 'FullLineageRanks']] = [Akkermansia, Akkermansia_rank]
    # df.loc[df.TaxID == 3115151, ['FullLineage', 'FullLineageRanks']] = [Akkermansia, Akkermansia_rank]
    # df.loc[df.TaxID == 3111450, ['FullLineage', 'FullLineageRanks']] = [Pseudarthrobacter, Pseudarthrobacter_rank]
    # df.loc[df.TaxID == 3117733, ['FullLineage', 'FullLineageRanks']] = [Halalkalicoccus, Halalkalicoccus_rank]
    # df.loc[df.TaxID == 3108483, ['FullLineage', 'FullLineageRanks']] = [Pseudarthrobacter, Pseudarthrobacter_rank]
    # df.loc[df.TaxID == 3090665, ['FullLineage', 'FullLineageRanks']] = [Pseudomonas, Pseudomonas_rank]
    # df.loc[df.TaxID == 3109594, ['FullLineage', 'FullLineageRanks']] = [Rhodococcus, Rhodococcus_rank]
    # df.loc[df.TaxID == 3110231, ['FullLineage', 'FullLineageRanks']] = [Thiobacillus, Thiobacillus_rank]
    # df.loc[df.TaxID == 3110227, ['FullLineage', 'FullLineageRanks']] = [Caproicibacter, Caproicibacter_rank]    
    # df.loc[df.TaxID == 3104265, ['FullLineage', 'FullLineageRanks']] = [Pseudomonas, Pseudomonas_rank]    
    # df.loc[df.TaxID == 3112258, ['FullLineage', 'FullLineageRanks']] = [Xanthomonas, Xanthomonas_rank]
    # df.loc[df.TaxID == 3117553, ['FullLineage', 'FullLineageRanks']] = [Bacillus, Bacillus_rank]
    # df.loc[df.TaxID == 3112561, ['FullLineage', 'FullLineageRanks']] = [Serratia, Serratia_rank]
    # df.loc[df.TaxID == 3112419, ['FullLineage', 'FullLineageRanks']] = [Exiguobacterium, Exiguobacterium_rank]

    print(df[df.FullLineage.isna()].head())
    print(df[df.FullLineageRanks.isna()].head())

    df['FullLineage'] = df['FullLineage'].fillna(zero_lineage)
    df['FullLineageRanks'] = df['FullLineageRanks'].fillna(zero_ranks)

    df.loc[:, 'FullLineage'] = df.loc[:, 'FullLineage'].str.split(';')
    df.loc[:, 'FullLineageRanks'] = df.loc[:, 'FullLineageRanks'].str.replace("no rank;", "root;").str.split(';')

    all_taxonomy_columns = ['superkingdom_centrifuge', 'clade_centrifuge', 'kingdom_centrifuge',
       'phylum_centrifuge', 'class_centrifuge', 'order_centrifuge', 'family_centrifuge',
       'genus_centrifuge', 'species_centrifuge', 'no rank_centrifuge', 'TaxID']

    # print(df[df.FullLineage.isna()].head())
    # print(df[df.FullLineageRanks.isna()].head())

    df.loc[:, 'dict_lineage'] = df.apply(lambda x: {rank + "_centrifuge":taxon for rank, taxon in zip(x.FullLineageRanks, x.FullLineage)} | {'TaxID':x.TaxID}, axis=1)

    records = df['dict_lineage'].to_list()

    taxo_df = pd.DataFrame(records)

    taxo_df.loc[:,all_taxonomy_columns[:-1]] = taxo_df.loc[:,all_taxonomy_columns[:-1]].fillna('Not Defined')

    return pl.from_pandas(taxo_df[all_taxonomy_columns])

#########################################################

def get_taxonomy_centrifuge(table, metadata, taxodir):
    """
    Count the number of misannotation in the table
    
    Parameters
    ----------
    path2table : str
        The path to the table
    metadata : pl.DataFrame
        The metadata of ICTV
    taxodir : str
        The path to the taxonomy directory
    
    Returns
    -------
    pd.DataFrame
        The misannotation count
    """

    bacteria_in_db = ["AE006468", "CP015418", "CP000031", "CP000830", "CP001312", "CP001357", "BX897699"]
    
    table = table.filter(~pl.col('Species ID').is_in(bacteria_in_db))

    all_taxonomy_columns = ['Virus GENBANK accession', 'Superkingdom', 'Realm', 'Kingdom',
       'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Virus name(s)']

    metadata = metadata.select(all_taxonomy_columns)

    table = table.join(metadata, left_on='Species ID', right_on='Virus GENBANK accession', how="left")    

    taxo_df = pytaxonkit.lineage(table['Taxonomy ID'].unique(), data_dir=taxodir)

    taxo_df = transform_taxo_centrifuge(taxo_df)

    taxo_df = taxo_df.fill_null('Not Defined')

    taxo_df = taxo_df.with_columns(
        pl.col('TaxID').cast(pl.Utf8)
    )

    table = table.with_columns(
        pl.col('Taxonomy ID').cast(pl.Utf8)
    )

    table = table.join(taxo_df, left_on='Taxonomy ID', right_on="TaxID", how="left")

    # table = table.filter(pl.col('superkingdom_kraken2') == 'Viruses')

    table = table.with_columns(
                    pl.col('File Name').str.replace(r'.centrifuge.classify', '')
                )

    all_taxonomy_columns = ['superkingdom_centrifuge', 'clade_centrifuge', 'kingdom_centrifuge',
        'phylum_centrifuge', 'class_centrifuge', 'order_centrifuge', 'family_centrifuge',
        'genus_centrifuge', 'species_centrifuge', 'no rank_centrifuge']

    # table = table.with_columns(
    #     [pl.col(column).fill_null('Not Defined').alias(column) for column in all_taxonomy_columns],
    # )

    # table = table.with_columns(
    #     [pl.col(column).fill_nan('Not Defined').alias(column) for column in all_taxonomy_columns]
    # )

    return table

#########################################################

def transform_taxo_kraken2(df):
    """
    Transform the taxo DataFrame
    
    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to be transformed
    
    Returns
    -------
    pl.DataFrame
        The transformed DataFrame
    """
    
    zero_ranks = 'superkingdom;clade;kingdom;phylum;class;order;family;subfamily;genus;species;no rank'
    zero_lineage = 'Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined'

    df.loc[df.TaxID == 0, ['FullLineage', 'FullLineageRanks']] = [zero_lineage, zero_ranks]
    df.loc[df.TaxID == 1, ['FullLineage', 'FullLineageRanks']] = [zero_lineage, zero_ranks]

    print(df[df.FullLineage.isna()].head())
    print(df[df.FullLineageRanks.isna()].head())

    df['FullLineage'] = df['FullLineage'].fillna(zero_lineage)
    df['FullLineageRanks'] = df['FullLineageRanks'].fillna(zero_ranks)

    df.loc[:, 'FullLineage'] = df.loc[:, 'FullLineage'].str.split(';')
    df.loc[:, 'FullLineageRanks'] = df.loc[:, 'FullLineageRanks'].str.split(';')

    all_taxonomy_columns = ['superkingdom_kraken2', 'clade_kraken2', 'kingdom_kraken2',
       'phylum_kraken2', 'class_kraken2', 'order_kraken2', 'family_kraken2',
       'genus_kraken2', 'species_kraken2', 'no rank_kraken2', 'TaxID']

    df.loc[:, 'dict_lineage'] = df.apply(lambda x: {rank + "_kraken2":taxon for rank, taxon in zip(x.FullLineageRanks, x.FullLineage)} | {'TaxID':x.TaxID}, axis=1)

    records = df['dict_lineage'].to_list()

    taxo_df = pd.DataFrame(records)

    taxo_df.loc[:,all_taxonomy_columns[:-1]] = taxo_df.loc[:,all_taxonomy_columns[:-1]].fillna('Not Defined')

    return pl.from_pandas(taxo_df[all_taxonomy_columns])

#########################################################

def get_taxonomy_kraken2(table, metadata, taxodir):
    """
    Count the number of misannotation in the table
    
    Parameters
    ----------
    path2table : str
        The path to the table
    metadata : pl.DataFrame
        The metadata of ICTV
    taxodir : str
        The path to the taxonomy directory
    
    Returns
    -------
    pd.DataFrame
        The misannotation count
    """
    # table = pl.read_csv(path2table, separator='\t')

    bacteria_in_db = ["AE006468", "CP015418", "CP000031", "CP000830", "CP001312", "CP001357", "BX897699"]
    
    table = table.filter(~pl.col('Species ID').is_in(bacteria_in_db))

    all_taxonomy_columns = ['Virus GENBANK accession', 'Superkingdom', 'Realm', 'Kingdom',
       'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Virus name(s)']

    metadata = metadata.select(all_taxonomy_columns)

    table = table.join(metadata, left_on='Species ID', right_on='Virus GENBANK accession', how="left")    

    taxo_df = pytaxonkit.lineage(table['Taxonomy ID'].unique(), data_dir=taxodir)

    taxo_df = transform_taxo_kraken2(taxo_df)

    taxo_df = taxo_df.fill_null('Not Defined')

    taxo_df = taxo_df.with_columns(
        pl.col('TaxID').cast(pl.Utf8)
    )

    table = table.with_columns(
        pl.col('Taxonomy ID').cast(pl.Utf8)
    )

    table = table.join(taxo_df, left_on='Taxonomy ID', right_on="TaxID", how="left")

    table = table.with_columns(
                    pl.col('File Name').str.replace(r'.kraken2.classify', '')
                )

    all_taxonomy_columns = ['superkingdom_kraken2', 'clade_kraken2', 'kingdom_kraken2',
       'phylum_kraken2', 'class_kraken2', 'order_kraken2', 'family_kraken2',
       'genus_kraken2', 'species_kraken2', 'no rank_kraken2']

    # table = table.with_columns(
    #     [pl.col(column).fill_null('Not Defined').alias(column) for column in all_taxonomy_columns],
    # )

    # table = table.with_columns(
    #     [pl.col(column).fill_nan('Not Defined').alias(column) for column in all_taxonomy_columns]
    # )


    return table

#########################################################

def transform_taxo_metabuli(df):
    """
    Transform the taxo DataFrame
    
    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to be transformed
    
    Returns
    -------
    pl.DataFrame
        The transformed DataFrame
    """
    
    zero_ranks = 'superkingdom;clade;kingdom;phylum;class;order;family;subfamily;genus;species;no rank'
    zero_lineage = 'Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined'
    # sopEPhi = 'Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Peduoviridae;Felsduovirus;Felsduovirus SopEphi'
    # sopEPhi_rank = 'superkingdom;clade;kingdom;phylum;class;family;genus;species'
    # medusavirus = 'Viruses;Varidnaviria;Bamfordvirae;Nucleocytoviricota;Megaviricetes;Mamonoviridae;Medusavirus;Medusavirus medusae'
    # medusavirus_rank = 'superkingdom;clade;kingdom;phylum;class;family;genus;species'
    # mischivirus = 'Viruses;Riboviria;Orthornavirae;Pisuviricota;Pisoniviricetes;Picornavirales;Picornaviridae;Caphthovirinae;Mischivirus;Mischivirus E'
    # mischivirus_rank = 'superkingdom;clade;kingdom;phylum;class;order;family;subfamily;genus;species'
    # nudivirus = 'Viruses;Naldaviricetes;Lefavirales;Nudiviridae;Betanudivirus;Betanudivirus hezeae'
    # nudivirus_rank = 'superkingdom;clade;order;family;genus;species'
    # sicinivirus = 'Viruses;Riboviria;Orthornavirae;Pisuviricota;Pisoniviricetes;Picornavirales;Picornaviridae;Kodimesavirinae;Sicinivirus;Sicinivirus A'
    # sicinivirus_rank = 'superkingdom;clade;kingdom;phylum;class;order;family;subfamily;genus;species'
    # Maribacter = 'Bacteria;Bacteroidota;Flavobacteriia;Flavobacteriales;Flavobacteriaceae;Maribacter;unclassified Maribacter'
    # Maribacter_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Dyella = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Lysobacterales;Rhodanobacteraceae;Dyella;unclassified Dyella'
    # Dyella_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Limnochorda = 'Bacteria;Bacillota;Limnochordia;Limnochordales;Limnochordaceae;Limnochorda;unclassified Limnochorda'
    # Limnochorda_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Lentzea = 'Bacteria;Actinomycetota;Actinobacteria;Actinomycetes;Pseudonocardiales;Pseudonocardiaceae;Lentzea;unclassified Lentzea'
    # Lentzea_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Vagococcus = 'Bacteria;Bacillota;Bacilli;Lactobacillales;Enterococcaceae;Vagococcus;unclassified Vagococcus'
    # Vagococcus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Chryseobacterium = 'Bacteria;Bacteroidota;Flavobacteriia;Flavobacteriales;Weeksellaceae;Chryseobacterium;unclassified Chryseobacterium'
    # Chryseobacterium_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Sedimentibacter = 'Bacteria;Bacillota;Tissierellia;Sedimentibacteriales;Sedimentibacter;unclassified Sedimentibacter'
    # Sedimentibacter_rank = 'superkingdom;phylum;class;genus;species'
    # Paenibacillus_silvisoli = 'Bacteria;Bacillota;Bacilli;Bacillales;Paenibacillaceae;Paenibacillus;Paenibacillus silvisoli'
    # Paenibacillus_silvisoli_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Haloarcula = 'Archaea;Euryarchaeota;Halobacteria;Halobacteriales;Halobacteriaceae;Haloarcula;unclassified Haloarcula'
    # Haloarcula_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Lelliottia = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Lelliottia;unclassified Lelliottia'
    # Lelliottia_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Pseudomonas = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;unclassified Pseudomonas'
    # Pseudomonas_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Klebsiella = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Klebsiella;unclassified Klebsiella'
    # Klebsiella_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Pseudoalteromonas = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Alteromonadales;Pseudoalteromonadaceae;Pseudoalteromonas;unclassified Pseudoalteromonas'
    # Pseudoalteromonas_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Sinanaerobacter = 'Bacteria;Bacillota;Clostridia;Eubacteriales;Sinanaerobacter;unclassified Sinanaerobacter'
    # Sinanaerobacter_rank = 'superkingdom;phylum;class;order;genus;species'
    # Streptomyces = 'Bacteria;Actinomycetota;Actinomycetes;Kitasatosporales;Streptomycetaceae;Streptomyces;unclassified Streptomyces'
    # Streptomyces_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Thermococcus = 'Archaea;Euryarchaeota;Thermococci;Thermococcales;Thermococcaceae;Thermococcus;unclassified Thermococcus'
    # Thermococcus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Bacteroides = 'Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;unclassified Bacteroides'
    # Bacteroides_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Qipengyuania = 'Bacteria;Pseudomonadota;Alphaproteobacteria;Sphingomonadales;Erythrobacteraceae;Qipengyuania;unclassified Qipengyuania'
    # Qipengyuania_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Nguyenibacter = 'Bacteria;Pseudomonadota;Alphaproteobacteria;Rhodospirillales;Acetobacteraceae;Nguyenibacter;unclassified Nguyenibacter'
    # Nguyenibacter_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Grimontia = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Vibrionales;Vibrionaceae;Grimontia;unclassified Grimontia'
    # Grimontia_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Novosphingobium = 'Bacteria;Pseudomonadota;Alphaproteobacteria;Sphingomonadales;Sphingomonadaceae;Novosphingobium;unclassified Novosphingobium'
    # Novosphingobium_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Akkermansia = 'Bacteria;Verrucomicrobiota;Verrucomicrobiia;Verrucomicrobiales;Akkermansiaceae;Akkermansia;unclassified Akkermansia'
    # Akkermansia_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Pseudarthrobacter = 'Bacteria;Actinomycetota;Actinomycetes;Micrococcales;Micrococcaceae;Pseudarthrobacter;unclassified Pseudarthrobacter'
    # Pseudarthrobacter_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Halalkalicoccus = 'Archaea;Euryarchaeota;Halobacteria;Halobacteriales;Halococcaceae;Halalkalicoccus;unclassified Halalkalicoccus'
    # Halalkalicoccus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Mycoplasma = 'Bacteria;Mycoplasmota;Mollicutes;Mycoplasmatales;Mycoplasmataceae;Mycoplasma;unclassified Mycoplasma'
    # Mycoplasma_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Rhodococcus = 'Bacteria;Actinomycetota;Actinomycetes;Myobacteriales;Nocardiaceae;Rhodococcus;unclassified Rhodococcus'
    # Rhodococcus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Thiobacillus = 'Bacteria;Pseudomonadota;Betaproteobacteria;Nitrosomonadales;Thiobacillaceae;Thiobacillus;unclassified Thiobacillus'
    # Thiobacillus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Caproicibacter = 'Bacteria;Bacillota;Clostridia;Eubacteriales;Acutalibacteraceae;Caproicibacter;unclassified Caproicibacter'
    # Caproicibacter_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Xanthomonas = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Lysobacterales;Lysobacteraceae;Xanthomonas;unclassified Xanthomonas'
    # Xanthomonas_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Bacillus = 'Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;unclassified Bacillus'
    # Bacillus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Serratia = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Yersiniaceae;Serratia;unclassified Serratia'
    # Serratia_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Exiguobacterium = 'Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Exiguobacterium;unclassified Exiguobacterium'
    # Exiguobacterium_rank = 'superkingdom;phylum;class;order;family;genus;species'

    df.loc[df.TaxID == 0, ['FullLineage', 'FullLineageRanks']] = [zero_lineage, zero_ranks]
    df.loc[df.TaxID == 1, ['FullLineage', 'FullLineageRanks']] = [zero_lineage, zero_ranks]
    # df.loc[df.TaxID == 3118720, ['FullLineage', 'FullLineageRanks']] = [sopEPhi, sopEPhi_rank]
    # df.loc[df.TaxID == 3114988, ['FullLineage', 'FullLineageRanks']] = [medusavirus, medusavirus_rank]
    # df.loc[df.TaxID == 2870396, ['FullLineage', 'FullLineageRanks']] = [mischivirus, mischivirus_rank]
    # df.loc[df.TaxID == 3116536, ['FullLineage', 'FullLineageRanks']] = [nudivirus, nudivirus_rank]
    # df.loc[df.TaxID == 2848078, ['FullLineage', 'FullLineageRanks']] = [sicinivirus, sicinivirus_rank]
    # df.loc[df.TaxID == 3053613, ['FullLineage', 'FullLineageRanks']] = [Maribacter, Maribacter_rank]
    # df.loc[df.TaxID == 3069105, ['FullLineage', 'FullLineageRanks']] = [Dyella, Dyella_rank]
    # df.loc[df.TaxID == 3109564, ['FullLineage', 'FullLineageRanks']] = [Limnochorda, Limnochorda_rank]
    # df.loc[df.TaxID == 3108822, ['FullLineage', 'FullLineageRanks']] = [Lentzea, Lentzea_rank]
    # df.loc[df.TaxID == 3109030, ['FullLineage', 'FullLineageRanks']] = [Vagococcus, Vagococcus_rank]
    # df.loc[df.TaxID == 3116594, ['FullLineage', 'FullLineageRanks']] = [Chryseobacterium, Chryseobacterium_rank]
    # df.loc[df.TaxID == 3109366, ['FullLineage', 'FullLineageRanks']] = [Sedimentibacter, Sedimentibacter_rank]
    # df.loc[df.TaxID == 3110539, ['FullLineage', 'FullLineageRanks']] = [Paenibacillus_silvisoli, Paenibacillus_silvisoli_rank]
    # df.loc[df.TaxID == 3109565, ['FullLineage', 'FullLineageRanks']] = [Limnochorda, Limnochorda_rank]
    # df.loc[df.TaxID == 3111776, ['FullLineage', 'FullLineageRanks']] = [Haloarcula, Haloarcula_rank]
    # df.loc[df.TaxID == 3110110, ['FullLineage', 'FullLineageRanks']] = [Lelliottia, Lelliottia_rank]
    # df.loc[df.TaxID == 3110772, ['FullLineage', 'FullLineageRanks']] = [Pseudomonas, Pseudomonas_rank]
    # df.loc[df.TaxID == 3111629, ['FullLineage', 'FullLineageRanks']] = [Klebsiella, Klebsiella_rank]
    # df.loc[df.TaxID == 3110111, ['FullLineage', 'FullLineageRanks']] = [Pseudomonas, Pseudomonas_rank]
    # df.loc[df.TaxID == 3112573, ['FullLineage', 'FullLineageRanks']] = [Pseudoalteromonas, Pseudoalteromonas_rank]
    # df.loc[df.TaxID == 3111540, ['FullLineage', 'FullLineageRanks']] = [Sinanaerobacter, Sinanaerobacter_rank]
    # df.loc[df.TaxID == 3111774, ['FullLineage', 'FullLineageRanks']] = [Streptomyces, Streptomyces_rank]
    # df.loc[df.TaxID == 3111325, ['FullLineage', 'FullLineageRanks']] = [Thermococcus, Thermococcus_rank]
    # df.loc[df.TaxID == 3117552, ['FullLineage', 'FullLineageRanks']] = [Bacteroides, Bacteroides_rank]
    # df.loc[df.TaxID == 3113984, ['FullLineage', 'FullLineageRanks']] = [Qipengyuania, Qipengyuania_rank]
    # df.loc[df.TaxID == 3049350, ['FullLineage', 'FullLineageRanks']] = [Nguyenibacter, Nguyenibacter_rank]
    # df.loc[df.TaxID == 3111011, ['FullLineage', 'FullLineageRanks']] = [Grimontia, Grimontia_rank]
    # df.loc[df.TaxID == 3109595, ['FullLineage', 'FullLineageRanks']] = [Novosphingobium, Novosphingobium_rank]
    # df.loc[df.TaxID == 3115152, ['FullLineage', 'FullLineageRanks']] = [Akkermansia, Akkermansia_rank]
    # df.loc[df.TaxID == 3115151, ['FullLineage', 'FullLineageRanks']] = [Akkermansia, Akkermansia_rank]
    # df.loc[df.TaxID == 3111450, ['FullLineage', 'FullLineageRanks']] = [Pseudarthrobacter, Pseudarthrobacter_rank]
    # df.loc[df.TaxID == 3117733, ['FullLineage', 'FullLineageRanks']] = [Halalkalicoccus, Halalkalicoccus_rank]
    # df.loc[df.TaxID == 3108483, ['FullLineage', 'FullLineageRanks']] = [Pseudarthrobacter, Pseudarthrobacter_rank]
    # df.loc[df.TaxID == 3090665, ['FullLineage', 'FullLineageRanks']] = [Pseudomonas, Pseudomonas_rank]
    # df.loc[df.TaxID == 3109594, ['FullLineage', 'FullLineageRanks']] = [Rhodococcus, Rhodococcus_rank]
    # df.loc[df.TaxID == 3110231, ['FullLineage', 'FullLineageRanks']] = [Thiobacillus, Thiobacillus_rank]
    # df.loc[df.TaxID == 3110227, ['FullLineage', 'FullLineageRanks']] = [Caproicibacter, Caproicibacter_rank]    
    # df.loc[df.TaxID == 3104265, ['FullLineage', 'FullLineageRanks']] = [Pseudomonas, Pseudomonas_rank]    
    # df.loc[df.TaxID == 3112258, ['FullLineage', 'FullLineageRanks']] = [Xanthomonas, Xanthomonas_rank]
    # df.loc[df.TaxID == 3117553, ['FullLineage', 'FullLineageRanks']] = [Bacillus, Bacillus_rank]
    # df.loc[df.TaxID == 3112561, ['FullLineage', 'FullLineageRanks']] = [Serratia, Serratia_rank]
    # df.loc[df.TaxID == 3112419, ['FullLineage', 'FullLineageRanks']] = [Exiguobacterium, Exiguobacterium_rank]

    print(df[df.FullLineage.isna()].head())
    print(df[df.FullLineageRanks.isna()].head())

    df['FullLineage'] = df['FullLineage'].fillna(zero_lineage)
    df['FullLineageRanks'] = df['FullLineageRanks'].fillna(zero_ranks)

    df.loc[:, 'FullLineage'] = df.loc[:, 'FullLineage'].str.split(';')
    df.loc[:, 'FullLineageRanks'] = df.loc[:, 'FullLineageRanks'].str.replace("no rank;", "root;").str.split(';')

    all_taxonomy_columns = ['superkingdom_metabuli', 'clade_metabuli', 'kingdom_metabuli',
       'phylum_metabuli', 'class_metabuli', 'order_metabuli', 'family_metabuli',
       'genus_metabuli', 'species_metabuli', 'no rank_metabuli', 'TaxID']

    # print(df[df.FullLineage.isna()].head())
    # print(df[df.FullLineageRanks.isna()].head())

    df.loc[:, 'dict_lineage'] = df.apply(lambda x: {rank + "_metabuli":taxon for rank, taxon in zip(x.FullLineageRanks, x.FullLineage)} | {'TaxID':x.TaxID}, axis=1)

    records = df['dict_lineage'].to_list()

    taxo_df = pd.DataFrame(records)

    taxo_df.loc[:,all_taxonomy_columns[:-1]] = taxo_df.loc[:,all_taxonomy_columns[:-1]].fillna('Not Defined')

    return pl.from_pandas(taxo_df[all_taxonomy_columns])

#########################################################

def get_taxonomy_metabuli(table, metadata, taxodir):
    """
    Count the number of misannotation in the table
    
    Parameters
    ----------
    path2table : str
        The path to the table
    metadata : pl.DataFrame
        The metadata of ICTV
    taxodir : str
        The path to the taxonomy directory
    
    Returns
    -------
    pd.DataFrame
        The misannotation count
    """

    bacteria_in_db = ["AE006468", "CP015418", "CP000031", "CP000830", "CP001312", "CP001357", "BX897699"]
    
    table = table.filter(~pl.col('Species ID').is_in(bacteria_in_db))

    all_taxonomy_columns = ['Virus GENBANK accession', 'Superkingdom', 'Realm', 'Kingdom',
       'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Virus name(s)']

    metadata = metadata.select(all_taxonomy_columns)

    table = table.join(metadata, left_on='Species ID', right_on='Virus GENBANK accession', how="left")    

    taxo_df = pytaxonkit.lineage(table['Taxonomy ID'].unique(), data_dir=taxodir)

    taxo_df = transform_taxo_metabuli(taxo_df)

    taxo_df = taxo_df.fill_null('Not Defined')

    taxo_df = taxo_df.with_columns(
        pl.col('TaxID').cast(pl.Utf8)
    )

    table = table.with_columns(
        pl.col('Taxonomy ID').cast(pl.Utf8)
    )

    table = table.join(taxo_df, left_on='Taxonomy ID', right_on="TaxID", how="left")

    # table = table.filter(pl.col('superkingdom_kraken2') == 'Viruses')

    table = table.with_columns(
                    pl.col('File Name').str.replace(r'.metabuli.classify', '')
                )

    all_taxonomy_columns = ['superkingdom_metabuli', 'clade_metabuli', 'kingdom_metabuli',
        'phylum_metabuli', 'class_metabuli', 'order_metabuli', 'family_metabuli',
        'genus_metabuli', 'species_metabuli', 'no rank_metabuli']

    # table = table.with_columns(
    #     [pl.col(column).fill_null('Not Defined').alias(column) for column in all_taxonomy_columns],
    # )

    # table = table.with_columns(
    #     [pl.col(column).fill_nan('Not Defined').alias(column) for column in all_taxonomy_columns]
    # )

    return table

#########################################################

def transform_taxo_krakenuniq(df):
    """
    Transform the taxo DataFrame
    
    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to be transformed
    
    Returns
    -------
    pl.DataFrame
        The transformed DataFrame
    """
    
    zero_ranks = 'superkingdom;clade;kingdom;phylum;class;order;family;subfamily;genus;species;no rank'
    zero_lineage = 'Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined;Not Defined'
    # sopEPhi = 'Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Peduoviridae;Felsduovirus;Felsduovirus SopEphi'
    # sopEPhi_rank = 'superkingdom;clade;kingdom;phylum;class;family;genus;species'
    # medusavirus = 'Viruses;Varidnaviria;Bamfordvirae;Nucleocytoviricota;Megaviricetes;Mamonoviridae;Medusavirus;Medusavirus medusae'
    # medusavirus_rank = 'superkingdom;clade;kingdom;phylum;class;family;genus;species'
    # mischivirus = 'Viruses;Riboviria;Orthornavirae;Pisuviricota;Pisoniviricetes;Picornavirales;Picornaviridae;Caphthovirinae;Mischivirus;Mischivirus E'
    # mischivirus_rank = 'superkingdom;clade;kingdom;phylum;class;order;family;subfamily;genus;species'
    # nudivirus = 'Viruses;Naldaviricetes;Lefavirales;Nudiviridae;Betanudivirus;Betanudivirus hezeae'
    # nudivirus_rank = 'superkingdom;clade;order;family;genus;species'
    # sicinivirus = 'Viruses;Riboviria;Orthornavirae;Pisuviricota;Pisoniviricetes;Picornavirales;Picornaviridae;Kodimesavirinae;Sicinivirus;Sicinivirus A'
    # sicinivirus_rank = 'superkingdom;clade;kingdom;phylum;class;order;family;subfamily;genus;species'
    # Maribacter = 'Bacteria;Bacteroidota;Flavobacteriia;Flavobacteriales;Flavobacteriaceae;Maribacter;unclassified Maribacter'
    # Maribacter_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Dyella = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Lysobacterales;Rhodanobacteraceae;Dyella;unclassified Dyella'
    # Dyella_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Limnochorda = 'Bacteria;Bacillota;Limnochordia;Limnochordales;Limnochordaceae;Limnochorda;unclassified Limnochorda'
    # Limnochorda_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Lentzea = 'Bacteria;Actinomycetota;Actinobacteria;Actinomycetes;Pseudonocardiales;Pseudonocardiaceae;Lentzea;unclassified Lentzea'
    # Lentzea_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Vagococcus = 'Bacteria;Bacillota;Bacilli;Lactobacillales;Enterococcaceae;Vagococcus;unclassified Vagococcus'
    # Vagococcus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Chryseobacterium = 'Bacteria;Bacteroidota;Flavobacteriia;Flavobacteriales;Weeksellaceae;Chryseobacterium;unclassified Chryseobacterium'
    # Chryseobacterium_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Sedimentibacter = 'Bacteria;Bacillota;Tissierellia;Sedimentibacteriales;Sedimentibacter;unclassified Sedimentibacter'
    # Sedimentibacter_rank = 'superkingdom;phylum;class;genus;species'
    # Paenibacillus_silvisoli = 'Bacteria;Bacillota;Bacilli;Bacillales;Paenibacillaceae;Paenibacillus;Paenibacillus silvisoli'
    # Paenibacillus_silvisoli_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Haloarcula = 'Archaea;Euryarchaeota;Halobacteria;Halobacteriales;Halobacteriaceae;Haloarcula;unclassified Haloarcula'
    # Haloarcula_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Lelliottia = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Lelliottia;unclassified Lelliottia'
    # Lelliottia_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Pseudomonas = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;unclassified Pseudomonas'
    # Pseudomonas_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Klebsiella = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Klebsiella;unclassified Klebsiella'
    # Klebsiella_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Pseudoalteromonas = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Alteromonadales;Pseudoalteromonadaceae;Pseudoalteromonas;unclassified Pseudoalteromonas'
    # Pseudoalteromonas_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Sinanaerobacter = 'Bacteria;Bacillota;Clostridia;Eubacteriales;Sinanaerobacter;unclassified Sinanaerobacter'
    # Sinanaerobacter_rank = 'superkingdom;phylum;class;order;genus;species'
    # Streptomyces = 'Bacteria;Actinomycetota;Actinomycetes;Kitasatosporales;Streptomycetaceae;Streptomyces;unclassified Streptomyces'
    # Streptomyces_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Thermococcus = 'Archaea;Euryarchaeota;Thermococci;Thermococcales;Thermococcaceae;Thermococcus;unclassified Thermococcus'
    # Thermococcus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Bacteroides = 'Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;unclassified Bacteroides'
    # Bacteroides_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Qipengyuania = 'Bacteria;Pseudomonadota;Alphaproteobacteria;Sphingomonadales;Erythrobacteraceae;Qipengyuania;unclassified Qipengyuania'
    # Qipengyuania_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Nguyenibacter = 'Bacteria;Pseudomonadota;Alphaproteobacteria;Rhodospirillales;Acetobacteraceae;Nguyenibacter;unclassified Nguyenibacter'
    # Nguyenibacter_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Grimontia = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Vibrionales;Vibrionaceae;Grimontia;unclassified Grimontia'
    # Grimontia_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Novosphingobium = 'Bacteria;Pseudomonadota;Alphaproteobacteria;Sphingomonadales;Sphingomonadaceae;Novosphingobium;unclassified Novosphingobium'
    # Novosphingobium_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Akkermansia = 'Bacteria;Verrucomicrobiota;Verrucomicrobiia;Verrucomicrobiales;Akkermansiaceae;Akkermansia;unclassified Akkermansia'
    # Akkermansia_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Pseudarthrobacter = 'Bacteria;Actinomycetota;Actinomycetes;Micrococcales;Micrococcaceae;Pseudarthrobacter;unclassified Pseudarthrobacter'
    # Pseudarthrobacter_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Halalkalicoccus = 'Archaea;Euryarchaeota;Halobacteria;Halobacteriales;Halococcaceae;Halalkalicoccus;unclassified Halalkalicoccus'
    # Halalkalicoccus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Mycoplasma = 'Bacteria;Mycoplasmota;Mollicutes;Mycoplasmatales;Mycoplasmataceae;Mycoplasma;unclassified Mycoplasma'
    # Mycoplasma_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Rhodococcus = 'Bacteria;Actinomycetota;Actinomycetes;Myobacteriales;Nocardiaceae;Rhodococcus;unclassified Rhodococcus'
    # Rhodococcus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Thiobacillus = 'Bacteria;Pseudomonadota;Betaproteobacteria;Nitrosomonadales;Thiobacillaceae;Thiobacillus;unclassified Thiobacillus'
    # Thiobacillus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Caproicibacter = 'Bacteria;Bacillota;Clostridia;Eubacteriales;Acutalibacteraceae;Caproicibacter;unclassified Caproicibacter'
    # Caproicibacter_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Xanthomonas = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Lysobacterales;Lysobacteraceae;Xanthomonas;unclassified Xanthomonas'
    # Xanthomonas_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Bacillus = 'Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;unclassified Bacillus'
    # Bacillus_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Serratia = 'Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Yersiniaceae;Serratia;unclassified Serratia'
    # Serratia_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Exiguobacterium = 'Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Exiguobacterium;unclassified Exiguobacterium'
    # Exiguobacterium_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Haloferax_volcanii = 'Archaea;Euryarchaeota;Halobacteria;Halobacteriales;Haloferacaceae;Haloferax;Haloferax volcanii'
    # Haloferax_volcanii_rank = 'superkingdom;phylum;class;order;family;genus;species'
    # Botryosphaeria_dothidea_virus = 'Viruses;Riboviria;Polyploviridae;Polymycovirus;Botryosphaeria dothidea virus'
    # Botryosphaeria_dothidea_virus_rank = 'superkingdom;clade;family;genus;species'

    df.loc[df.TaxID == 0, ['FullLineage', 'FullLineageRanks']] = [zero_lineage, zero_ranks]
    df.loc[df.TaxID == 1, ['FullLineage', 'FullLineageRanks']] = [zero_lineage, zero_ranks]
    # df.loc[df.TaxID == 3118720, ['FullLineage', 'FullLineageRanks']] = [sopEPhi, sopEPhi_rank]
    # df.loc[df.TaxID == 3114988, ['FullLineage', 'FullLineageRanks']] = [medusavirus, medusavirus_rank]
    # df.loc[df.TaxID == 2870396, ['FullLineage', 'FullLineageRanks']] = [mischivirus, mischivirus_rank]
    # df.loc[df.TaxID == 3116536, ['FullLineage', 'FullLineageRanks']] = [nudivirus, nudivirus_rank]
    # df.loc[df.TaxID == 2848078, ['FullLineage', 'FullLineageRanks']] = [sicinivirus, sicinivirus_rank]
    # df.loc[df.TaxID == 3053613, ['FullLineage', 'FullLineageRanks']] = [Maribacter, Maribacter_rank]
    # df.loc[df.TaxID == 3069105, ['FullLineage', 'FullLineageRanks']] = [Dyella, Dyella_rank]
    # df.loc[df.TaxID == 3109564, ['FullLineage', 'FullLineageRanks']] = [Limnochorda, Limnochorda_rank]
    # df.loc[df.TaxID == 3108822, ['FullLineage', 'FullLineageRanks']] = [Lentzea, Lentzea_rank]
    # df.loc[df.TaxID == 3109030, ['FullLineage', 'FullLineageRanks']] = [Vagococcus, Vagococcus_rank]
    # df.loc[df.TaxID == 3116594, ['FullLineage', 'FullLineageRanks']] = [Chryseobacterium, Chryseobacterium_rank]
    # df.loc[df.TaxID == 3109366, ['FullLineage', 'FullLineageRanks']] = [Sedimentibacter, Sedimentibacter_rank]
    # df.loc[df.TaxID == 3110539, ['FullLineage', 'FullLineageRanks']] = [Paenibacillus_silvisoli, Paenibacillus_silvisoli_rank]
    # df.loc[df.TaxID == 3109565, ['FullLineage', 'FullLineageRanks']] = [Limnochorda, Limnochorda_rank]
    # df.loc[df.TaxID == 3111776, ['FullLineage', 'FullLineageRanks']] = [Haloarcula, Haloarcula_rank]
    # df.loc[df.TaxID == 3110110, ['FullLineage', 'FullLineageRanks']] = [Lelliottia, Lelliottia_rank]
    # df.loc[df.TaxID == 3110772, ['FullLineage', 'FullLineageRanks']] = [Pseudomonas, Pseudomonas_rank]
    # df.loc[df.TaxID == 3111629, ['FullLineage', 'FullLineageRanks']] = [Klebsiella, Klebsiella_rank]
    # df.loc[df.TaxID == 3110111, ['FullLineage', 'FullLineageRanks']] = [Pseudomonas, Pseudomonas_rank]
    # df.loc[df.TaxID == 3112573, ['FullLineage', 'FullLineageRanks']] = [Pseudoalteromonas, Pseudoalteromonas_rank]
    # df.loc[df.TaxID == 3111540, ['FullLineage', 'FullLineageRanks']] = [Sinanaerobacter, Sinanaerobacter_rank]
    # df.loc[df.TaxID == 3111774, ['FullLineage', 'FullLineageRanks']] = [Streptomyces, Streptomyces_rank]
    # df.loc[df.TaxID == 3111325, ['FullLineage', 'FullLineageRanks']] = [Thermococcus, Thermococcus_rank]
    # df.loc[df.TaxID == 3117552, ['FullLineage', 'FullLineageRanks']] = [Bacteroides, Bacteroides_rank]
    # df.loc[df.TaxID == 3113984, ['FullLineage', 'FullLineageRanks']] = [Qipengyuania, Qipengyuania_rank]
    # df.loc[df.TaxID == 3049350, ['FullLineage', 'FullLineageRanks']] = [Nguyenibacter, Nguyenibacter_rank]
    # df.loc[df.TaxID == 3111011, ['FullLineage', 'FullLineageRanks']] = [Grimontia, Grimontia_rank]
    # df.loc[df.TaxID == 3109595, ['FullLineage', 'FullLineageRanks']] = [Novosphingobium, Novosphingobium_rank]
    # df.loc[df.TaxID == 3115152, ['FullLineage', 'FullLineageRanks']] = [Akkermansia, Akkermansia_rank]
    # df.loc[df.TaxID == 3115151, ['FullLineage', 'FullLineageRanks']] = [Akkermansia, Akkermansia_rank]
    # df.loc[df.TaxID == 3111450, ['FullLineage', 'FullLineageRanks']] = [Pseudarthrobacter, Pseudarthrobacter_rank]
    # df.loc[df.TaxID == 3117733, ['FullLineage', 'FullLineageRanks']] = [Halalkalicoccus, Halalkalicoccus_rank]
    # df.loc[df.TaxID == 3108483, ['FullLineage', 'FullLineageRanks']] = [Pseudarthrobacter, Pseudarthrobacter_rank]
    # df.loc[df.TaxID == 3090665, ['FullLineage', 'FullLineageRanks']] = [Pseudomonas, Pseudomonas_rank]
    # df.loc[df.TaxID == 3109594, ['FullLineage', 'FullLineageRanks']] = [Rhodococcus, Rhodococcus_rank]
    # df.loc[df.TaxID == 3110231, ['FullLineage', 'FullLineageRanks']] = [Thiobacillus, Thiobacillus_rank]
    # df.loc[df.TaxID == 3110227, ['FullLineage', 'FullLineageRanks']] = [Caproicibacter, Caproicibacter_rank]    
    # df.loc[df.TaxID == 3104265, ['FullLineage', 'FullLineageRanks']] = [Pseudomonas, Pseudomonas_rank]    
    # df.loc[df.TaxID == 3112258, ['FullLineage', 'FullLineageRanks']] = [Xanthomonas, Xanthomonas_rank]
    # df.loc[df.TaxID == 3117553, ['FullLineage', 'FullLineageRanks']] = [Bacillus, Bacillus_rank]
    # df.loc[df.TaxID == 3112561, ['FullLineage', 'FullLineageRanks']] = [Serratia, Serratia_rank]
    # df.loc[df.TaxID == 3112419, ['FullLineage', 'FullLineageRanks']] = [Exiguobacterium, Exiguobacterium_rank]
    # df.loc[df.TaxID == 114529, ['FullLineage', 'FullLineageRanks']] = [Haloferax_volcanii, Haloferax_volcanii_rank]
    # df.loc[df.TaxID == 1516075, ['FullLineage', 'FullLineageRanks']] = [Botryosphaeria_dothidea_virus, Botryosphaeria_dothidea_virus_rank]

    print(df[df.FullLineage.isna()].head())
    print(df[df.FullLineageRanks.isna()].head())

    df['FullLineage'] = df['FullLineage'].fillna(zero_lineage)
    df['FullLineageRanks'] = df['FullLineageRanks'].fillna(zero_ranks)

    df.loc[:, 'FullLineage'] = df.loc[:, 'FullLineage'].str.split(';')
    df.loc[:, 'FullLineageRanks'] = df.loc[:, 'FullLineageRanks'].str.replace("no rank;", "root;").str.split(';')

    all_taxonomy_columns = ['superkingdom_krakenuniq', 'clade_krakenuniq', 'kingdom_krakenuniq',
       'phylum_krakenuniq', 'class_krakenuniq', 'order_krakenuniq', 'family_krakenuniq',
       'genus_krakenuniq', 'species_krakenuniq', 'no rank_krakenuniq', 'TaxID']


    df.loc[:, 'dict_lineage'] = df.apply(lambda x: {rank + "_krakenuniq":taxon for rank, taxon in zip(x.FullLineageRanks, x.FullLineage)} | {'TaxID':x.TaxID}, axis=1)

    records = df['dict_lineage'].to_list()

    taxo_df = pd.DataFrame(records)

    taxo_df.loc[:,all_taxonomy_columns[:-1]] = taxo_df.loc[:,all_taxonomy_columns[:-1]].fillna('Not Defined')

    return pl.from_pandas(taxo_df[all_taxonomy_columns])

#########################################################

def get_taxonomy_krakenuniq(table, metadata, taxodir):
    """
    Count the number of misannotation in the table
    
    Parameters
    ----------
    path2table : str
        The path to the table
    metadata : pl.DataFrame
        The metadata of ICTV
    taxodir : str
        The path to the taxonomy directory
    
    Returns
    -------
    pd.DataFrame
        The misannotation count
    """

    bacteria_in_db = ["AE006468", "CP015418", "CP000031", "CP000830", "CP001312", "CP001357", "BX897699"]
    
    table = table.filter(~pl.col('Species ID').is_in(bacteria_in_db))

    all_taxonomy_columns = ['Virus GENBANK accession', 'Superkingdom', 'Realm', 'Kingdom',
       'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Virus name(s)']

    metadata = metadata.select(all_taxonomy_columns)

    table = table.join(metadata, left_on='Species ID', right_on='Virus GENBANK accession', how="left")    

    taxo_df = pytaxonkit.lineage(table['Taxonomy ID'].unique(), data_dir=taxodir)

    taxo_df = transform_taxo_krakenuniq(taxo_df)

    taxo_df = taxo_df.fill_null('Not Defined')

    taxo_df = taxo_df.with_columns(
        pl.col('TaxID').cast(pl.Utf8)
    )

    table = table.with_columns(
        pl.col('Taxonomy ID').cast(pl.Utf8)
    )

    table = table.join(taxo_df, left_on='Taxonomy ID', right_on="TaxID", how="left")

    # table = table.filter(pl.col('superkingdom_kraken2') == 'Viruses')

    table = table.with_columns(
                    pl.col('File Name').str.replace(r'.krakenuniq.classify', '')
                )

    all_taxonomy_columns = ['superkingdom_krakenuniq', 'clade_krakenuniq', 'kingdom_krakenuniq',
        'phylum_krakenuniq', 'class_krakenuniq', 'order_krakenuniq', 'family_krakenuniq',
        'genus_krakenuniq', 'species_krakenuniq', 'no rank_krakenuniq']

    # table = table.with_columns(
    #     [pl.col(column).fill_null('Not Defined').alias(column) for column in all_taxonomy_columns],
    # )

    # table = table.with_columns(
    #     [pl.col(column).fill_nan('Not Defined').alias(column) for column in all_taxonomy_columns]
    # )

    return table



#########################################################
# 3. Gobal variables
#########################################################

all_bacteria_ICTV = [
    "AE006468",
    "AF049230",
    "BD143114",
    "BD269513",
    "BX897699",
    "CP000031",
    "CP000830",
    "CP001312",
    "CP001357",
    "CP006891",
    "CP014526",
    "CP015418",
    "CP019275",
    "CP023680",
    "CP023686",
    "CP038625",
    "CP052639",
    "JAEILC010000038",
    "KK213166",
    "M18706",
    "QUVN01000024",
    "U68072",
    "U96748",
]
provirus_not_well_defined = [
    "AY319521",
    "EF462197",
    "EF462198",
    "EF710638",
    "EF710642",
    "FJ184280",
    "FJ188381",
    "HG424323",
    "HM208303",
    "J02013",
    "JQ347801",
    "K02712",
    "KF147927",
    "KF183314",
    "KF183315",
    "KP972568",
    "KX232515",
    "KX452695",
    "MK075003",
    "V01201",
]
unverified = [
    "DQ188954",
    "HM543472",
    "JQ407224",
    "KC008572",
    "KC626021",
    "KF302037",
    "KF360047",
    "KF938901",
    "KJ641726",
    "KM233624",
    "KM389459",
    "KM982402",
    "KP843857",
    "KR862307",
    "KU343148",
    "KU343149",
    "KU343150",
    "KU343151",
    "KU343152",
    "KU343153",
    "KU343154",
    "KU343155",
    "KU343156",
    "KU343160",
    "KU343161",
    "KU343162",
    "KU343163",
    "KU343164",
    "KU343165",
    "KU343169",
    "KU343170",
    "KU343171",
    "KU672593",
    "KU752557",
    "KX098515",
    "KX228196",
    "KX228197",
    "KX228198",
    "KX363561",
    "KX452696",
    "KX452698",
    "KX656670",
    "KX656671",
    "KX989546",
    "KY450753",
    "KY487839",
    "KY608967",
    "KY742649",
    "MG459218",
    "MG551742",
    "MG599035",
    "MH791395",
    "MH791402",
    "MH791405",
    "MH791410",
    "MH791412",
    "MH918795",
    "MH925094",
    "MH992121",
    "MK033136",
    "MK050014",
    "MK415316",
    "MK415317",
    "MK474470",
    "MK780203",
    "MN545971",
    "MN871450",
    "MN871491",
    "MN871495",
    "MN871498",
    "MN928506",
    "MT360681",
    "MT360682",
    "MW325771",
    "MW685514",
    "MW685515",
]
big_unwanted_list = all_bacteria_ICTV + provirus_not_well_defined + unverified

#########################################################

print("Start the script")

print("Read the metadata")
metadata = pl.read_csv(snakemake.input.metadata, separator='\t')
metadata = metadata.filter(~pl.col('Virus GENBANK accession').is_in(big_unwanted_list))

# Remove if not bacteria
metadata = metadata.with_columns(
    pl.lit("bacteria").alias("Kingdom"),
    pl.lit("bacteria").alias("Virus name(s)"),
    pl.lit("bacteria").alias("Realm"),
)

# If virus
# metadata = metadata.with_columns(
#     pl.lit("Viruses").alias("Superkingdom")
# )

print(metadata.columns)

print("Read the output file")
if snakemake.params.taxo_profiler == "centrifuge":
    output_centrifuge = Path(snakemake.input.tsv)

    print("Parse the output file - Centrifuge")
    tmp = parse_output_centrifuge_files(output_centrifuge, big_unwanted_list)

    print("Get the taxonomy - Centrifuge")
    tmp = get_taxonomy_centrifuge(
        tmp,
        metadata,
        '/mnt/archgen/microbiome_coprolite/aVirus/03-data/refdbs/metabuli_db/taxonomy'
        # "/mnt/archgen/microbiome_coprolite/aVirus/03-data/refdbs/ICTV/centrifuge_ICTV/taxonomy"
    )
elif snakemake.params.taxo_profiler == "kraken2":
    output_kraken2 = Path(snakemake.input.tsv)

    print("Parse the output file - Kraken2")
    tmp = parse_output_kraken2_files(output_kraken2, big_unwanted_list) 

    print("Get the taxonomy - Kraken2")
    tmp = get_taxonomy_kraken2(
        tmp,
        metadata,
        '/mnt/archgen/microbiome_coprolite/aVirus/03-data/refdbs/metabuli_db/taxonomy'
        # "/mnt/archgen/microbiome_coprolite/aVirus/03-data/refdbs/ICTV/kraken2_ICTV/taxonomy"
    )
elif snakemake.params.taxo_profiler == "metabuli":
    output_metabuli = Path(snakemake.input.tsv)

    print("Parse the output file - Metabuli")
    tmp = parse_output_metabuli_files(output_metabuli, big_unwanted_list)

    print("Get the taxonomy - Metabuli")
    tmp = get_taxonomy_metabuli(
        tmp,
        metadata,
        '/mnt/archgen/microbiome_coprolite/aVirus/03-data/refdbs/metabuli_db/taxonomy'
        # "/mnt/archgen/microbiome_coprolite/aVirus/03-data/refdbs/ICTV/metabuli_ICTV/taxonomy"
    )
elif snakemake.params.taxo_profiler == "krakenuniq":
    output_krakenuniq = Path(snakemake.input.tsv)

    print("Parse the output file - Krakenuniq")
    # Same format as Kraken2
    tmp = parse_output_krakenuniq_files(output_krakenuniq, big_unwanted_list)

    print("Get the taxonomy - Krakenuniq")
    tmp = get_taxonomy_krakenuniq(
        tmp,
        metadata,
        '/mnt/archgen/microbiome_coprolite/aVirus/03-data/refdbs/metabuli_db/taxonomy'
        # "/mnt/archgen/microbiome_coprolite/aVirus/03-data/refdbs/ICTV/krakenuniq_ICTV/taxonomy"
    )

if snakemake.params.taxo_profiler == "centrifuge":
    rank = [
        ("superkingdom_centrifuge", "Superkingdom"),
        ("clade_centrifuge", "Realm"),
        ("kingdom_centrifuge", "Kingdom"),
        ("phylum_centrifuge", "Phylum"),
        ("class_centrifuge", "Class"),
        # ("order_centrifuge", "Order"), 
        # ("family_centrifuge", "Family"), 
        ("genus_centrifuge", "Genus"), 
        ("species_centrifuge", "Species"), 
        # ("no rank_centrifuge", "Virus name(s)"), 
        ]
elif snakemake.params.taxo_profiler == "kraken2":
    rank = [
        ("superkingdom_kraken2", "Superkingdom"),
        ("clade_kraken2", "Realm"),
        ("kingdom_kraken2", "Kingdom"),
        ("phylum_kraken2", "Phylum"),
        ("class_kraken2", "Class"),
        # ("order_kraken2", "Order"), 
        # ("family_kraken2", "Family"), 
        ("genus_kraken2", "Genus"), 
        ("species_kraken2", "Species"), 
        # ("no rank_kraken2", "Virus name(s)"), 
        ]
elif snakemake.params.taxo_profiler == "metabuli":
    rank = [
        ("superkingdom_metabuli", "Superkingdom"),
        ("clade_metabuli", "Realm"),
        ("kingdom_metabuli", "Kingdom"),
        ("phylum_metabuli", "Phylum"),
        ("class_metabuli", "Class"),
        # ("order_metabuli", "Order"), 
        # ("family_metabuli", "Family"), 
        ("genus_metabuli", "Genus"), 
        ("species_metabuli", "Species"), 
        # ("no rank_metabuli", "Virus name(s)"), 
        ]
elif snakemake.params.taxo_profiler == "krakenuniq":
    rank = [
        ("superkingdom_krakenuniq", "Superkingdom"),
        ("clade_krakenuniq", "Realm"),
        ("kingdom_krakenuniq", "Kingdom"),
        ("phylum_krakenuniq", "Phylum"),
        ("class_krakenuniq", "Class"),
        # ("order_krakenuniq", "Order"), 
        # ("family_krakenuniq", "Family"), 
        ("genus_krakenuniq", "Genus"), 
        ("species_krakenuniq", "Species"), 
        # ("no rank_krakenuniq", "Virus name(s)"), 
        ]

print("Transform the taxonomy")
# If virus please change it back to normal this 
tmp = tmp.with_columns(
    # pl.lit("Viruses").alias("Superkingdom"),
    pl.lit("No rank equal").alias("equal_rank"),
    pl.lit(False).alias("equal_higher"),
)


for profiler_rank, ictv_rank in rank:
    if ictv_rank == "Species":
        tmp = tmp.with_columns(
            pl.col(profiler_rank).eq(pl.col(ictv_rank)).alias("equal_species"),
        ).with_columns(
            pl.col('equal_species').replace_strict(
                    old=True,
                    new="Species",
                    default=pl.col('equal_rank'),
            ).alias("equal_rank"),
        )
    else :
        tmp = tmp.with_columns(
            pl.col(profiler_rank).eq(pl.col(ictv_rank)).alias("equal_higher_tmp"),
        ).with_columns( # We used a tmp column to avoid the conflict with the next column as we go from the higher rank to the lower rank
            pl.col('equal_higher_tmp').replace_strict( 
                    old=True,
                    new=True,
                    default=pl.col('equal_higher'),
            ).alias("equal_higher"),
        ).with_columns(
            pl.col('equal_higher_tmp').replace_strict(
                    old=True,
                    new=ictv_rank,
                    default=pl.col('equal_rank'),
            ).alias("equal_rank"),
        )

tmp = tmp.drop('equal_higher_tmp')

if snakemake.params.taxo_profiler == "centrifuge":
    columns_profiler_superkingdom = "superkingdom_centrifuge"
elif snakemake.params.taxo_profiler == "kraken2":
    columns_profiler_superkingdom = "superkingdom_kraken2"
elif snakemake.params.taxo_profiler == "metabuli":
    columns_profiler_superkingdom = "superkingdom_metabuli"
elif snakemake.params.taxo_profiler == "krakenuniq":
    columns_profiler_superkingdom = "superkingdom_krakenuniq"

tmp = tmp.with_columns(
    pl.col(columns_profiler_superkingdom).eq(pl.col('Superkingdom')).alias("equal_superkingdom"),
)

tmp = tmp.with_columns(
    pl.col(columns_profiler_superkingdom).replace_strict(
            old='Not Defined',
            new='Unclassified',
            default=pl.col('equal_higher'),
    ).alias("equal_higher"),
    pl.col(columns_profiler_superkingdom).replace_strict(
            old='Not Defined',
            new='Unclassified',
            default=pl.col('equal_species'),
    ).alias("equal_species"),
    pl.col(columns_profiler_superkingdom).replace_strict(
            old='Not Defined',
            new='Unclassified',
            default=pl.col('equal_superkingdom'),
    ).alias("equal_superkingdom"),
) 


tmp = tmp.filter(
    ~pl.col('Species ID').is_in(big_unwanted_list)
)


tmp = tmp.with_columns(
    pl.col('equal_higher').fill_null(
        pl.col('equal_superkingdom')
    ),
)
# tmp = tmp.with_columns(
#     pl.col('equal_higher').fill_nan(
#         pl.col('equal_superkingdom')
#     ),
# )


# tmp.to_pandas().to_csv(snakemake.output.taxonomy, index=False, sep='\t')
tmp.write_parquet(snakemake.output.taxonomy)


    


