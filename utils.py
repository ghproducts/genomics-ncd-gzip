
import re
import os
import pandas as pd
import numpy as np


class TaxaLevel:
    def __init__(self, taxa_dir:str):
        self.names = None
        self.nodes = None
        self.taxa_dir = taxa_dir

    def parse_files(self):
        """
        Parses the files in the given directory and places the names
        and nodes in the TaxaLevel object.

        Args:
            taxa_dir (str): The directory containing the taxa files.
        
        Thie looks in the taxa_dir for the files names.dmp and nodes.dmp
        """
        files = os.listdir(self.taxa_dir)
        if 'nodes.dmp' not in files:
            raise FileNotFoundError("nodes.dmp not found in the specified directory.")
        elif 'names.dmp' not in files:
            print("names.dmp not found in the specified directory. This may cause issues with taxa parsing in the future.")
 
        node_filename = os.path.join(self.taxa_dir, 'nodes.dmp')
        df = pd.read_csv(node_filename, sep='|', header=None, engine='c',
                         usecols=[0, 1, 2], names=['taxid', 'parent', 'rank'],
                         dtype={'taxid': int, 'parent': int, 'rank': str})
        df['rank'] = df['rank'].str.strip()
        self.nodes = df.set_index('taxid')

        if 'names.dmp' in files:
            names_filename = os.path.join(self.taxa_dir, 'names.dmp')
            df = pd.read_csv(names_filename, sep='|', header=None, engine='c',
                            usecols=[0, 1], names=['taxid', 'name'],
                            dtype={'taxid': int, 'name': str})
            df['name'] = df['name'].str.strip()
            self.names = df.set_index('taxid')['name'].to_dict()


    def link_taxid(self, taxid:int, level:str, return_name=False):
        """
        returns the higher-level taxid of a genome at a specified level.

        Args:
            taxid (int): The taxid of the genome.
            level (str): The level to which the taxid should be linked.
            return_name (bool): If True, returns the name of taxa as wekk as the taxid.

        returns:
            int or tuple: The taxid of the higher-level taxa at the specified level.
                          If return_name is True, returns a tuple of (taxid, name).
        """

        if self.nodes is None:
            raise ValueError("TaxaLevel object has not been initialized with nodes data.")
        if (self.names is None) and return_name:
            raise ValueError("TaxaLevel object has not been initialized with names data. make sure to parse the names.dmp file if requesting taxa names")
        if level not in self.nodes['rank'].values:
            raise ValueError(f"Level '{level}' not found in taxa ranks.")
        
        
        TAXONOMIC_HIERARCHY = [
            'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family',
            'genus', 'species', 'subspecies', 'no rank'
        ]

        current_rank = self.nodes.loc[taxid]['rank']
        if (
            current_rank in TAXONOMIC_HIERARCHY and
            level in TAXONOMIC_HIERARCHY and
            TAXONOMIC_HIERARCHY.index(level) >= TAXONOMIC_HIERARCHY.index(current_rank)
        ):
            raise ValueError(f"Requested level '{level}' must be higher than current level '{current_rank}'.")


        # begin parsing the taxa tree
        current = taxid
        visited = set()

        while current != 1 and current not in visited:
            visited.add(current)
            try:
                row = self.nodes.loc[current]
            except KeyError:
                # If the taxid is not found in the nodes, return None
                print(f"Taxid {current} not found in nodes.")
                return None

            if row['rank'] == level:
                if return_name:
                    name = self.names.get(current, None)
                    return (current, name)
                else:
                    return current
                
            current = int(row['parent'])

        return None


valid_seq = re.compile(r'^[ACGT]+$')  # compiled regex is very fast

def sample_valid_fragments(seq, n, fragment_size):
    fragments = []
    max_start = len(seq) - fragment_size
    attempts = 0
    while len(fragments) < n and attempts < n * 5:
        start = np.random.randint(0, max_start)
        fragment = seq[start:start + fragment_size]
        if valid_seq.fullmatch(fragment):
            fragments.append(fragment)
        attempts += 1
    return fragments

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))