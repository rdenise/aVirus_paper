#########################################################
# 1. Importing libraries
#########################################################

import polars as pl
import numpy as np

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
from seaborn import plotting_context

#########################################################
# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")


#########################################################
# 2. Functions
#########################################################

# def calculate_precision(column):
#     nb_reads = (column != 'Unclassified').sum()

#     if nb_reads == 0:
#         return np.nan

#     nb_correct = (column == 'true').sum()

#     return nb_correct / nb_reads

# #########################################################

def create_dataframe_mean(list_of_files, big_unwanted_list, read2size):
    tmp = [pl.read_parquet(file) for file in list_of_files]
    df = pl.concat(tmp, how="vertical_relaxed")

    # df = pl.read_parquet(
    #     list_of_files,
    # )

    df = df.filter(
        ~pl.col('Species ID').is_in(big_unwanted_list)
    )

    df = df.with_columns(
        pl.col('Read ID').str.split('-2').list[0],
    )

    df = df.with_columns(
        pl.col('Read ID').replace_strict(reads2size).alias("Read Length"),
    )

    df = df.with_columns(
        equal_genus = pl.col('equal_rank').replace_strict(
            {
                'Genus': 'true',
                'Species': 'true',
                'No rank equal': 'Unclassified',
            },
            default='false'
        )
    )

    df = df.filter(
        pl.col('equal_species') != 'Unclassified'
    )

    df_species = df.group_by('Read Length').agg(
        pl.col('equal_species').eq("true").sum() / pl.col('equal_species').count(),
    )

    df_genus = df.group_by('Read Length').agg(
            pl.col('equal_genus').eq("true").sum() / pl.col('equal_genus').count(),
    )

    df_higher = df.group_by('Read Length').agg(
        pl.col('equal_higher').eq("true").sum() / pl.col('equal_higher').count(),
    )

    df = df_species.join(df_higher, on='Read Length')
    df = df.join(df_genus, on='Read Length')

    df = df.rename(
        {
            "equal_species": "Species",
            "equal_higher": "Superkingdom",
            "equal_genus": "Genus",
        }
    )

    df = df.unpivot(
        on = ['Species', 'Superkingdom', 'Genus'],
        index = ['Read Length'],
        variable_name = 'Level',
        value_name = 'Precision',
    )

    df = df.sort(['Level', 'Read Length'], descending=[False, False])

    conditions, software, *leftover = Path(list_of_files[0]).stem.split('.')

    size, leftover = conditions.split('-lognormal-')
    damage, coverage = leftover.split('-deamination-')

    df = df.with_columns(
        pl.lit(damage.capitalize()).alias("Damage"),
        pl.lit(software.capitalize()).alias("Software"),
    )

    return df

#########################################################

def plot_precision_length(df):

    dict_damage_software = {
        'Centrifuge - Heavy': '#267DC1',
        'Centrifuge - Light': '#267DC1',
        'Centrifuge - None': '#267DC1',
        'Kraken2 - Heavy': '#FF8D2D',
        'Kraken2 - Light': '#FF8D2D',
        'Kraken2 - None': '#FF8D2D',
        'Metabuli - Heavy': '#14B92D',
        'Metabuli - Light': '#14B92D',
        'Metabuli - None': '#14B92D',
        'Krakenuniq - Heavy': '#FF2D2D',
        'Krakenuniq - Light': '#FF2D2D',
        'Krakenuniq - None': '#FF2D2D',
    }

    print('Plotting')

    # so.Plot.config.theme.update(
    #     axes_style("ticks",rc=custom_style),
    # )

    p = (
        so.Plot(df, x="Read Length", y="Precision", color="damage_software", linestyle="Level")
        .layout(size=(12, 9)) # width, height
        .add(so.Line(), so.Agg())
        .scale(
            color=dict_damage_software, 
            linestyle={'Species': '--', 'Superkingdom': '-', 'Genus': '-.'},
        )
        .limit(y=(0,1.01), x=(35, None))
        # .limit(y=(0,1.01), x=(35, 500))
        .theme(axes_style("whitegrid") | plotting_context("talk") | {"grid.linestyle": ":"})
    )

    print('Saving')

    p.save(snakemake.output.pdf, bbox_inches='tight')
    p.save(snakemake.output.png, bbox_inches='tight')

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

dict_files = {
    'kraken2': [],
    'centrifuge': [],
    'metabuli': [],
    'krakenuniq': [],
}

for taxo in all_taxo:
    conditions, software, *leftover = Path(taxo).stem.split('.')

    dict_files[software].append(taxo)

reads2size = pl.read_parquet(
    dict_files['kraken2']
)

reads2size = reads2size.with_columns(
    pl.col('Read ID').str.split('-2').list[0],
    size=pl.col('Read ID').str.split(':').list[-1].str.split('-2').list[0].cast(pl.Int32),
)

reads2size = dict(reads2size[['Read ID', 'size']].iter_rows())


df_mean = []

if Path(snakemake.params.tsv).exists():
    df_mean = pl.read_parquet(snakemake.params.tsv)
else:
    for software, files in dict_files.items():
        print(f'Doing: {software}')
        
        df_mean.append(
            create_dataframe_mean(files, big_unwanted_list,reads2size)
        )

    df_mean = pl.concat(df_mean)

    df_mean.write_parquet(snakemake.params.tsv)

df_mean = df_mean.with_columns(
    damage_software=pl.col('Software').str.to_titlecase() + ' - ' + pl.col('Damage').str.to_titlecase()
)

df_mean = df_mean.filter(
    pl.col('Read Length') <= 200
)

# p = plot_precision_length(df_mean[df_mean['Read Length'].isin(range(35, 501, 10))])
p = plot_precision_length(df_mean.filter(pl.col('Read Length').is_in(range(35, 200, 5))))