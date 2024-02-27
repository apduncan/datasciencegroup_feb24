"""Some quick function to import data from different profilers as abundance
tables."""

import glob
import os
from typing import Iterable
import pandas as pd


def read_metaphlan(f: os.PathLike, rank: str = "g") -> pd.DataFrame:
    """Read 

    Args:
        f (os.PathLike): _description_
        rank (str, optional): _description_. Defaults to "g".

    Returns:
        pd.DataFrame: MetaPhlAn format table at rank.
    """

    rank_letter: str = rank[0].lower()
    sample_name: str = os.path.basename(f)
    return (
        pd.read_csv(
            f,
            comment="#",
            names=["clade_name", "ncbi_tax_id", "relative_abundance",
                "additional_species"],
            index_col=0,
            sep="\t"
        )
        .filter(
            regex=rf'^.*{rank_letter}__[a-zA-Z0-9_]*$', 
            axis=0
        )
    )

def read_metaphlans(
        files: Iterable[os.PathLike],
        rank: str = "g") -> pd.DataFrame:
    """Combine multiple MetaPhlAn profiles to a single table at a given rank.

    Args:
        files (os.PathLike): Files for the sample profiles
        rank (str, optional): Rank to extract. Defaults to "g".

    Returns:
        pd.DataFrame: taxa x samples table of relative abundance at a specified
        taxonomic rank.
    """

    df: pd.DataFrame = pd.concat(
        (read_metaphlan(x)['relative_abundance'] for x in files),
        axis=1
    )
    df.columns = [os.path.basename(x) for x in files]
    return df

def read_kmcp(
        f: os.PathLike, 
        rank: str = "genus") -> pd.Series:
    df: pd.DataFrame = pd.read_csv(
        f,
        comment="@",
        names=["taxid", "rank", "taxpath", "taxpathsn", "percentage"],
        sep="\t"
    )
    # Convert taxpaths to use ; delimiters
    df['taxpathsn'] = df['taxpathsn'].apply(lambda x: x.replace('|', ';'))
    df = df.set_index('taxpathsn')
    df = df.loc[df['rank'] == rank]
    return df['percentage']

def read_kmcps(
        files: Iterable[os.PathLike],
        rank: str = "genus"
) -> pd.DataFrame:
    """Table of relative abundances based on KMCP cami.profile files"""
    df: pd.DataFrame = pd.concat(
        (read_kmcp(x, rank=rank) for x in files),
        axis=1
    )
    df.columns = [os.path.basename(x).split(".")[0] for x in files]
    return df

if __name__ == "__main__":
    # Dirty tests
    # mp: pd.DataFrame = read_metaphlan(
    #     "data/metaphlan/SRX5707173_R1.fastq.gz.tsv", rank="genus")
    # read_metaphlans(glob.glob("data/metaphlan/*.tsv"))

    kmcp: pd.Series = read_kmcp("data/kmcp/SRX5707173.cami.profile")
    kmcp_df: pd.DataFrame = read_kmcps(glob.glob("data/kmcp/*.cami.profile"))