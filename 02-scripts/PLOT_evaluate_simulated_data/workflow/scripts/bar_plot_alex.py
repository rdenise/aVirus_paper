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

def count_rank_assignation(taxo_path, profiler, big_unwanted_list):
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

    if profiler == "kraken2":
        rank = [
            ("no rank_kraken2", "Virus name(s)"), 
            ("species_kraken2", "Species"), 
            ("genus_kraken2", "Genus"), 
            ("class_kraken2", "Class"),
            ("phylum_kraken2", "Phylum"),
            ("kingdom_kraken2", "Kingdom"),
            ("clade_kraken2", "Realm"),
            ("superkingdom_kraken2", "Superkingdom"),
            ]
    elif profiler == "centrifuge":
        rank = [
            ("no rank_centrifuge", "Virus name(s)"), 
            ("species_centrifuge", "Species"), 
            ("genus_centrifuge", "Genus"), 
            ("class_centrifuge", "Class"),
            ("phylum_centrifuge", "Phylum"),
            ("kingdom_centrifuge", "Kingdom"),
            ("clade_centrifuge", "Realm"),
            ("superkingdom_centrifuge", "Superkingdom"),
            ]
    elif profiler == "metabuli":
        rank = [
            ("no rank_metabuli", "Virus name(s)"), 
            ("species_metabuli", "Species"), 
            ("genus_metabuli", "Genus"), 
            ("class_metabuli", "Class"),
            ("phylum_metabuli", "Phylum"),
            ("kingdom_metabuli", "Kingdom"),
            ("clade_metabuli", "Realm"),
            ("superkingdom_metabuli", "Superkingdom"),
            ]
    elif profiler == "krakenuniq":
        rank = [
            ("no rank_krakenuniq", "Virus name(s)"), 
            ("species_krakenuniq", "Species"), 
            ("genus_krakenuniq", "Genus"), 
            ("class_krakenuniq", "Class"),
            ("phylum_krakenuniq", "Phylum"),
            ("kingdom_krakenuniq", "Kingdom"),
            ("clade_krakenuniq", "Realm"),
            ("superkingdom_krakenuniq", "Superkingdom"),
            ]

    tmp = [pl.read_parquet(file) for file in taxo_path]
    # taxo_df = pl.read_parquet(taxo_path)
    taxo_df = pl.concat(tmp, how="vertical_relaxed")

    

    taxo_df = taxo_df.filter(
            ~pl.col('Species ID').is_in(big_unwanted_list)
        )

    taxo_df = taxo_df.with_columns(
        pl.col("File Name").str.split('-').list[0].str.to_titlecase().alias("Size"),
        pl.col("File Name").str.split('-').list[2].str.to_titlecase().alias("Deamination"),
        pl.col("File Name").str.split('-').list[4].replace_strict(num2cov).alias("Coverage"),
    )

    taxo_df = taxo_df.with_columns(
        Conditions = pl.col("Size") + " - " + pl.col("Deamination") + " deamination - " + pl.col("Coverage")
    )

    # Get the total reads per condition
    total_reads = taxo_df.group_by('Conditions').len(name='Reads')

    list_rank = []

    for profiler_rank, ictv_rank in rank:
        # Get the rank assignation
        count_df = taxo_df.with_columns(
            pl.col(profiler_rank).eq(pl.col(ictv_rank)).alias("x == y"),
        )

        count_df = count_df.with_columns(
            pl.col(profiler_rank).replace_strict(
                old="Not Defined",
                new="Not Defined",
                default=pl.col("x == y"),
            ).alias("x == y"),
        ).group_by(['Conditions', 'x == y']).len(name=ictv_rank)

        count_df = count_df.with_columns(pl.col("x == y").fill_null("Not Defined"))

        count_df = count_df.join(total_reads, on='Conditions', how='left')

        count_df = count_df.with_columns(
            (pl.col(ictv_rank) / pl.col("Reads") * 100).round(2)
        )

        count_df = count_df.with_columns(
            tmp = pl.col("Conditions").str.split(' - ')
        ).with_columns(
            Size=pl.col("tmp").list[0],
            Deamination=pl.col("tmp").list[1].str.replace(" deamination", ""),
            Coverage=pl.col("tmp").list[2],
            Coverage_sort=pl.col("tmp").list[2].str.replace("x", "").cast(pl.Float32),
        ).drop(['Reads', "tmp"])

        list_rank.append(count_df)
    
    full_count_melt = pl.concat(list_rank, how='align').sort(['Deamination','Coverage_sort']).drop(['Coverage_sort'])

    full_count_melt = full_count_melt.fill_null(0)

    full_count_melt = full_count_melt.unpivot(
        index=["Size", "Deamination", "Coverage", "Conditions", "x == y"], 
        on=cs.numeric(), 
        variable_name="Rank", 
        value_name="Percentage")

    return full_count_melt

#########################################################

def plot_true_false(df):
    """
    Plot a bar chart showing the percentage of true and false values in a DataFrame.

    Parameters:
    -----------
    df : pandas.DataFrame
        The DataFrame containing the data to be plotted.

    Returns:
    --------
    seaborn.objects.Plot 
        The plotted bar chart.

    """
    
    so.Plot.config.theme.update(
        axes_style("ticks",rc=custom_style | {"axes.spines.bottom": False}),
    )

    df = df.with_columns(
        pl.col("x == y").replace_strict(
            {
                "Not Defined":"aaNot Defined", 
                "true":"true", 
                "false":"false"
            }
        )
    ).sort(["x == y", "Rank", "Profiler"], descending=True)

    df = df.with_columns(
        pl.col("x == y").replace_strict(
            {
                "aaNot Defined":"Not Defined",
                "true":"true",
                "false":"false"
            }
        ),
    )
    
    p = (
        so.Plot(df, x="Profiler", y="Percentage", color="x == y")
        .layout(size=(12, 9)) # width, height
        .facet("Rank", "Deamination",  order={"row": ["None", "Light", "Heavy"], "col": ["Virus name(s)", "Species", "Genus", "Class", "Superkingdom"]})
        .add(so.Bar(), so.Stack())
        .scale(color=so.Nominal(values=["#7aa457", "#a80321","#777777"],
                                order=["true", "false", "Not Defined"]),
            )
        .limit(y=(0, 100))
        .plot()
    )

    for ax in p._figure.axes :
        ax.xaxis.set_tick_params(rotation=90)

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
# 4. Main
#########################################################

all_taxo = sorted(snakemake.input.taxonomy)

count_melt = []

dict_files = {
    "kraken2": [],
    "centrifuge": [],
    "metabuli": [],
    "krakenuniq": [],
}

for file in all_taxo:
    profiler = Path(file).stem.split('.')[1]
    dict_files[profiler].append(file)

for profiler, files in dict_files.items():
    print(profiler)
    count_melt.append(count_rank_assignation(
                    files, 
                    profiler,
                    big_unwanted_list,
                )
    )

    count_melt[-1] = count_melt[-1].with_columns(
        Profiler = pl.lit(profiler),
    )

print("Concatenating")
count_melt = pl.concat(count_melt)

count_melt = count_melt.filter(pl.col('Rank').is_in(['Virus name(s)', 'Species', 'Genus', 'Class', 'Superkingdom']))

print("Plotting")
p = plot_true_false(count_melt)

p.save(snakemake.output.pdf, bbox_inches='tight')

plt.close()

# count_melt.to_pandas().to_csv(snakemake.output.tsv, sep='\t', index=False)
count_melt.write_parquet(snakemake.output.tsv)

