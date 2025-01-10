#########################################################
# 1. Importing libraries
#########################################################

import polars as pl
import polars.selectors as cs
import pandas as pd
from pathlib import Path
import sys

import matplotlib.pyplot as plt

custom_style = {"text.color": "#131516",
                "svg.fonttype": "none",
                # "font.family": "sans-serif",
                # "font.weight": "light",
                "axes.spines.right": False,
                "axes.spines.top": False,
                # "axes.spines.bottom": False,
                'xtick.bottom': False,
                }

import seaborn as sns ; sns.set_theme(style="ticks", rc=custom_style)  # for plot styling
from seaborn import axes_style
import seaborn.objects as so


#########################################################
# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")


#########################################################
# 2. Functions
#########################################################

def count_reads_annotation(taxonomy_df, rank, kraken_rank, size, info_df):
    # Count reads attributed to a specific species
    correct = taxonomy_df.group_by('Species ID').agg(
        pl.col(rank).eq(pl.col(kraken_rank)).alias("x == y").fill_null(False).sum()
                )


    # Count reads with correct annotation
    count_reads = taxonomy_df.group_by('Species ID').agg(pl.col("Species ID").count().alias("count"))

    # Create a new dataframe with the counts
    correct = correct.join(count_reads, on="Species ID", how="left")

    
    correct = correct.join(info_df.select(['Accession Number', 'minimizers']), right_on="Accession Number", left_on='Species ID', how="left")

    correct = correct.with_columns(
        percentage_true = pl.col("x == y") / pl.col("count") * 100
    )
    
    correct = correct.join(size_df, on="Species ID", how="left")

    return correct

#########################################################

def rank_reads_annotation(taxo_path, size_df, info_kraken, profiler, big_unwanted_list):
    """
    Count the rank assignation

    Parameters
    ----------
    taxo_path : str
        Path to the table with the taxonomy

    Returns
    -------
    pl.DataFrame
        The DataFrame with the rank assignation
    """

    num2cov = {'0':'0.1x', '1':'0.5x', '2':'1x', '3':'3x', '4':'5x', '5':'10x'}

    taxo_df =  pl.read_parquet(taxo_path)

    taxo_df = taxo_df.filter(
        ~pl.col('Species ID').is_in(big_unwanted_list)
    )

    taxo_df = taxo_df.with_columns(
        pl.lit("Viruses").alias("Superkingdom")
    )

    full_count_df = []  

    if profiler == "kraken2":
        rank = [
            ("no rank_kraken2", "Virus name(s)"), 
            ("species_kraken2", "Species"), 
            ("genus_kraken2", "Genus"), 
            ("class_kraken2", "Class"),
            # ("phylum_kraken2", "Phylum"),
            # ("kingdom_kraken2", "Kingdom"),
            # ("clade_kraken2", "Realm"),
            # ("superkingdom_kraken2", "Superkingdom"),
            ]
    elif profiler == "centrifuge":
        rank = [
            ("no rank_centrifuge", "Virus name(s)"), 
            ("species_centrifuge", "Species"), 
            ("genus_centrifuge", "Genus"), 
            ("class_centrifuge", "Class"),
            # ("phylum_centrifuge", "Phylum"),
            # ("kingdom_centrifuge", "Kingdom"),
            # ("clade_centrifuge", "Realm"),
            # ("superkingdom_centrifuge", "Superkingdom"),
            ]
    elif profiler == "metabuli":
        rank = [
            ("no rank_metabuli", "Virus name(s)"), 
            ("species_metabuli", "Species"), 
            ("genus_metabuli", "Genus"), 
            ("class_metabuli", "Class"),
            # ("phylum_metabuli", "Phylum"),
            # ("kingdom_metabuli", "Kingdom"),
            # ("clade_metabuli", "Realm"),
            # ("superkingdom_metabuli", "Superkingdom"),
            ]
    elif profiler == "krakenuniq":
        rank = [
            ("no rank_krakenuniq", "Virus name(s)"), 
            ("species_krakenuniq", "Species"), 
            ("genus_krakenuniq", "Genus"), 
            ("class_krakenuniq", "Class"),
            # ("phylum_krakenuniq", "Phylum"),
            # ("kingdom_krakenuniq", "Kingdom"),
            # ("clade_krakenuniq", "Realm"),
            # ("superkingdom_krakenuniq", "Superkingdom"),
            ]

    for profiler_rank, ictv_rank in rank:
        full_count_df.append(count_reads_annotation(
            taxo_df,
            ictv_rank,
            profiler_rank,
            size_df,
            info_kraken,
        ).with_columns(
            Rank = pl.lit(ictv_rank)
        )
        )

    full_count_df = pl.concat(full_count_df)

    conditions = Path(taxo_path).stem.split('.')[0]

    size, leftover = conditions.split('-lognormal-')
    damage, coverage = leftover.split('-deamination-')

    conditions = f"{size.capitalize()} - {damage.capitalize()} deamination - {num2cov[coverage]}"

    full_count_df = full_count_df.with_columns(
        Conditions = pl.lit(conditions),
        Coverage = pl.lit(num2cov[coverage])
        )

    return full_count_df

#########################################################

def plot_percentage_correct(correct, y_column="minimizers"):
   
    so.Plot.config.theme.update(
        axes_style("whitegrid"),
    )

    p = (
        so.Plot(correct, y=y_column, x="percentage_true")
        .layout(size=(12, 9)) # width, height
        .facet(col="Rank", row="Conditions")
        .add(so.Dots())
        .limit(x=(0, 100))
        .scale(y="log")
    )

    return p

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

long_none = [
    'Long - None deamination - 0.1x',
    'Long - None deamination - 0.5x',
    'Long - None deamination - 1x',
    'Long - None deamination - 3x',
    'Long - None deamination - 5x',
    'Long - None deamination - 10x',
]

long_light = [
    'Long - Light deamination - 0.1x',
    'Long - Light deamination - 0.5x',
    'Long - Light deamination - 1x',
    'Long - Light deamination - 3x',
    'Long - Light deamination - 5x',
    'Long - Light deamination - 10x',
]

long_heavy = [
    'Long - Heavy deamination - 0.1x',
    'Long - Heavy deamination - 0.5x',
    'Long - Heavy deamination - 1x',
    'Long - Heavy deamination - 3x',
    'Long - Heavy deamination - 5x',
    'Long - Heavy deamination - 10x',
 ]

medium_none = [
    'Medium - None deamination - 0.1x',
    'Medium - None deamination - 0.5x',
    'Medium - None deamination - 1x',
    'Medium - None deamination - 3x',
    'Medium - None deamination - 5x',
    'Medium - None deamination - 10x',
]

medium_light = [
    'Medium - Light deamination - 0.1x',
    'Medium - Light deamination - 0.5x',
    'Medium - Light deamination - 1x',
    'Medium - Light deamination - 3x',
    'Medium - Light deamination - 5x',
    'Medium - Light deamination - 10x',
]

medium_heavy = [
    'Medium - Heavy deamination - 0.1x',
    'Medium - Heavy deamination - 0.5x',
    'Medium - Heavy deamination - 1x',
    'Medium - Heavy deamination - 3x',
    'Medium - Heavy deamination - 5x',
    'Medium - Heavy deamination - 10x',
 ]

short_none = [
    'Short - None deamination - 0.1x',
    'Short - None deamination - 0.5x',
    'Short - None deamination - 1x',
    'Short - None deamination - 3x',
    'Short - None deamination - 5x',
    'Short - None deamination - 10x',
]
short_light = [
    'Short - Light deamination - 0.1x',
    'Short - Light deamination - 0.5x',
    'Short - Light deamination - 1x',
    'Short - Light deamination - 3x',
    'Short - Light deamination - 5x',
    'Short - Light deamination - 10x',
]
short_heavy = [
    'Short - Heavy deamination - 0.1x',
    'Short - Heavy deamination - 0.5x',
    'Short - Heavy deamination - 1x',
    'Short - Heavy deamination - 3x',
    'Short - Heavy deamination - 5x',
    'Short - Heavy deamination - 10x',
 ]

#########################################################

dict_size_deamination = {
    "Short - None": short_none,
    "Short - Light": short_light,
    "Short - Heavy": short_heavy,
    "Medium - None": medium_none,
    "Medium - Light": medium_light,
    "Medium - Heavy": medium_heavy,
    "Long - None": long_none,
    "Long - Light": long_light,
    "Long - Heavy": long_heavy,
}

#########################################################
# 4. Main
#########################################################

all_taxo = sorted(snakemake.input.taxonomy)

size_df = {"Species ID":[], "Size":[]}

with open("/mnt/archgen/microbiome_coprolite/aVirus/03-data/refdbs/ICTV/ICTV_database/20240209/ICTV_all_genomes.fna", "rt") as r_file:
    for line in r_file:
        if line.startswith(">"):
            species_id, size = line[1:].split()[:2]
            size = int(size[:-2])

            if species_id not in big_unwanted_list:
                size_df["Species ID"].append(species_id)
                size_df["Size"].append(size)

size_df = pl.DataFrame(size_df)

kraken_info =  pl.read_csv(snakemake.input.inspect, separator="\t")

kraken_info = kraken_info.join(
    pl.read_csv("/mnt/archgen/microbiome_coprolite/aVirus/03-data/refdbs/ICTV/ICTV_database/20240209/accession_taxonomic_table.csv"), 
    right_on="Taxonomic ID", left_on="taxonomy_id", how="left",
)

kraken_info = kraken_info.rename({"count":"minimizers"})

kraken_info = kraken_info.with_columns(
    pl.col("Accession Number").str.split(".").list[0]
)

kraken_info = kraken_info.filter(~pl.col('Accession Number').is_in(big_unwanted_list))

correct_melt = []

for file in all_taxo:
    correct_melt.append(
        rank_reads_annotation(
            file,
            size_df,
            kraken_info,
            snakemake.params.taxo_profiler,
            big_unwanted_list,
        )
    )
    
correct_melt = pl.concat(correct_melt)

p = plot_percentage_correct(correct_melt)

p.save(snakemake.output.pdf, bbox_inches='tight')
plt.close()

# correct_melt.to_pandas().to_csv(snakemake.output.tsv, sep='\t', index=False)
correct_melt.write_parquet(snakemake.output.tsv)