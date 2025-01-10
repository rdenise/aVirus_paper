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
                'pdf.fonttype': 42, # Editable text
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

def plot_sensitivity_VS_precision(df, name):

    # dict_color = {
    #     'Centrifuge - Heavy': '#267DC1',
    #     'Centrifuge - Light': '#74B1E0',
    #     'Centrifuge - None': '#CFE8FF',
    #     'Kraken2 - Heavy': '#FF8D2D',
    #     'Kraken2 - Light': '#F5AD71',
    #     'Kraken2 - None': '#F5D3B4',
    #     'Metabuli - Heavy': '#14B92D',
    #     'Metabuli - Light': '#63D773',
    #     'Metabuli - None': '#BBF2C4',
    #     'Krakenuniq - Heavy': '#FF2D2D',
    #     'Krakenuniq - Light': '#F57373',
    #     'Krakenuniq - None': '#F2C4C4',
    #     'Centrifuge - Short': '#267DC1', # Depending on the thing we want to highlight on the plot
    #     'Centrifuge - Medium': '#74B1E0',
    #     'Centrifuge - Long': '#CFE8FF',
    #     'Kraken2 - Short': '#FF8D2D',
    #     'Kraken2 - Medium': '#F5AD71',
    #     'Kraken2 - Long': '#F5D3B4',
    #     'Metabuli - Short': '#14B92D',
    #     'Metabuli - Medium': '#63D773',
    #     'Metabuli - Long': '#BBF2C4',
    #     'Krakenuniq - Short': '#FF2D2D',
    #     'Krakenuniq - Medium': '#F57373',
    #     'Krakenuniq - Long': '#F2C4C4',
    # }

    dict_color = {
        'Centrifuge - Species': '#267DC1',
        'Centrifuge - Genus': '#74B1E0',
        'Kraken2 - Species': '#FF8D2D',
        'Kraken2 - Genus': '#F5AD71',
        'Metabuli - Species': '#14B92D',
        'Metabuli - Genus': '#63D773',
        'Krakenuniq - Species': '#FF2D2D',
        'Krakenuniq - Genus': '#F57373',
    }

    print('Plotting')

    # so.Plot.config.theme.update(
    #     axes_style("ticks",rc=custom_style),
    # )

    p = (
        # so.Plot(df, x="Sensitivity", y="Precision", color=name, marker="Level")
        so.Plot(df, x="Sensitivity", y="Precision", color=name, marker="Size")
        .layout(size=(8, 8)) # width, height
        .add(so.Dots(pointsize=10, stroke=1.5, fillalpha=0.5))
        .scale(
            color=dict_color, 
            # marker={'Species': 'o', 'Species + Higher': '^', 'Species + Genus': 's'},
            marker={'Short': 'o', 'Medium': '^', 'Long': 's'},
        )
        .limit(x=(0,1.01), y=(0, 1.01))
        .theme(axes_style("whitegrid") | plotting_context("talk") | {"grid.linestyle": ":", "pdf.fonttype": 42})
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

# tmp = [pl.read_parquet(file) for file in sorted(snakemake.input.taxonomy)]
# all_taxo = pl.concat(tmp, how="vertical_relaxed")

all_taxo = sorted(snakemake.input.taxonomy)

df_mean = {
    'Sensitivity': [],
    'Precision': [],
    'Level': [],
    'Damage': [],
    'Size': [],
    'Coverage': [],
    'Software': [],
}

if Path(snakemake.params.tsv).exists():
    df_mean = pl.read_parquet(snakemake.params.tsv)
else:
    for taxo in all_taxo:
        print(f'Doing: {taxo}')
        df = pl.read_parquet(taxo)

        df = df.filter(
            ~pl.col('Species ID').is_in(big_unwanted_list)
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

        # Calculate sensitivity
        mean_sensitivity_species = df.group_by('Species ID').agg(
            pl.col('equal_species').eq("true").sum() / pl.col('equal_species').count(),
        )['equal_species'].mean()

        mean_sensitivity_genus = df.group_by('Species ID').agg(
            pl.col('equal_genus').eq("true").sum() / pl.col('equal_genus').count(),
        )['equal_genus'].mean()

        mean_sensitivity_higher = df.group_by('Species ID').agg(
            pl.col('equal_higher').eq("true").sum() / pl.col('equal_higher').count(),
        )['equal_higher'].mean()

        # Calculate precision
        df = df.filter(
            pl.col('equal_species') != 'Unclassified'
        )

        mean_precision_species = df.group_by('Species ID').agg(
            pl.col('equal_species').eq("true").sum() / pl.col('equal_species').count(),
        )['equal_species'].mean()

        mean_precision_genus = df.group_by('Species ID').agg(
            pl.col('equal_genus').eq("true").sum() / pl.col('equal_genus').count(),
        )['equal_genus'].mean()

        mean_precision_higher = df.group_by('Species ID').agg(
            pl.col('equal_higher').eq("true").sum() / pl.col('equal_higher').count(),
        )['equal_higher'].mean()

        conditions, software, *leftover = Path(taxo).stem.split('.')

        size, leftover = conditions.split('-lognormal-')
        damage, coverage = leftover.split('-deamination-')

        num2cov = {'0':'0.1x', '1':'0.5x', '2':'1x', '3':'3x', '4':'5x', '5':'10x'}

        df_mean['Sensitivity'] += [mean_sensitivity_species, mean_sensitivity_genus, mean_sensitivity_higher]
        df_mean['Precision'] += [mean_precision_species, mean_precision_genus, mean_precision_higher]
        df_mean['Level'] += ['Species', 'Species + Genus', 'Species + Higher']
        df_mean['Damage'] += [damage.capitalize()] * 3
        df_mean['Size'] += [size.capitalize()] * 3
        df_mean['Coverage'] += [num2cov[coverage]] * 3
        df_mean['Software'] += [software] * 3

    df_mean = pl.from_dict(df_mean)

    df_mean.write_parquet(snakemake.params.tsv)

if snakemake.params.focus == 'damage':
    df_mean = df_mean.with_columns(
        damage_software=pl.col('Software').str.to_titlecase() + ' - ' + pl.col('Damage').str.to_titlecase()
    )
    name_column = 'Damage Software'
elif snakemake.params.focus == 'size':
    df_mean = df_mean.with_columns(
        (pl.col('Software').str.to_titlecase() + ' - ' + pl.col('Size').str.to_titlecase()).alias('Software - Size'),
    )
    name_column = 'Software - Size'
elif snakemake.params.focus == 'level':
    df_mean = df_mean.with_columns(
        (pl.col('Software').str.to_titlecase() + ' - ' + pl.col('Level').str.split(' + ').list[-1]).alias('Software - Level'),
    )
    name_column = 'Software - Level'

    df_mean = df_mean.filter(
        pl.col('Level') != 'Species + Higher'
    )

p = plot_sensitivity_VS_precision(df_mean, name=name_column)