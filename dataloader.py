from dataclasses import dataclass
from typing import Literal, Optional
from Bio import SeqIO
from utils import TaxaLevel, sample_valid_fragments

import numpy as np
import os

TaxRank = Literal["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]

class Fasta_DataLoader:
    '''
    class to load data from a directory
    '''

    def __init__(self, path):
        self.path = path
        self.data = None # structured array 

    def load_data(self):
        '''
        function to load all fasta files from a directory
        Every fasta file within the directory will be taken and imported

        the fasta ID will be used for the genome label
        '''
        labels = []
        seqs = []

        for file in os.listdir(self.path):
            if file.endswith(('.fa', '.fasta', '.fna')):
                for record in SeqIO.parse(os.path.join(self.path, file), 'fasta'):
                    labels.append(record.id)
                    seqs.append(str(record.seq).upper())
        
        max_len = max(len(seq) for seq in seqs)
        dtype = np.dtype([
            ('label', 'U50'),
            ('seq', f'U{max_len}'),  # dynamic based on longest sequence
            ('clen', 'i4')
        ])
        self.data = np.empty(len(labels), dtype=dtype)
        self.data['label'] = labels
        self.data['seq'] = seqs
        self.data['clen'] = np.zeros(len(seqs))  # compressed length

    def subsample(self, 
                  mode: Literal["genome", "taxa"] = "genome",
                  fragment_size: Optional[int] = 126,
                  *,
                  S: Optional[int] = 1000,
                  rank: Optional[TaxRank] = 'species',
                  taxa_data: Optional[TaxaLevel] = None,
                  N: Optional[int] = 10):
        '''
        function to subsample genomes in a directory by taxa levels or genome size
        either select genome leve (default) or taxa level

        taxa-levels: 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'genome'
        subsampling by taxa results N samples from each taxa, down to species level. 
        N = samples per class in taxa
        
        if taxa-level is selected, lineage information must be provided. Also, the fasta
        header for each genome must contain the lineage information in the format:
        >taxID
        SEQUENCE...

        subsampling by genome level just splits each genome into fragments. The allocation of these
        fragments is porportional to the length of each genome
        N = L(training_genome)/L(fragment_size) * S
        S = consstant to reduce the porportion of samples from the genome

        args: 
            mode - genome or taxa
            fragment_size - size of the fragments to split the genomes into (default: 126)
            
        Genome level args:
            S - discussed earlier

        Taxa level args:
            rank - species, genus, family, order, class, phylum, kingdom 
            taxa_data - TaxaLevel object containing the taxa information
            N - discussed earlier
        '''
        if self.data is None:
            raise ValueError("Data not loaded. Please load data first.")

        dtype = np.dtype([
                    ('label', 'U50'),
                    ('seq', f'U{fragment_size}'),  # dynamic based on longest sequence
                    ('clen', 'i4')
                ])
        
        labels = []
        fragments = []

        if mode == "genome":
            # Subsample by genome level
            if S <= 0:
                raise ValueError("S must be a positive integer.")
            
            # Calculate number of samples per genome
            lengths = [len(seq) for seq in self.data['seq']]
            max_samples = [max(1, int(length / fragment_size)) for length in lengths]

            for idx, seq in enumerate(self.data['seq']):
                if len(seq) < fragment_size:
                    continue
                N_fragments = max_samples[idx] * S
                valid_fragments = sample_valid_fragments(seq, N_fragments, fragment_size)
                fragments.extend(valid_fragments)
                labels.extend([self.data['label'][idx]] * len(valid_fragments))

        elif mode == "taxa":
            #subsample by taxa level
            if taxa_data is None:
                raise ValueError("Taxa data must be provided for taxa-level subsampling.")
            
            taxa_at_rank = []
            for taxid in self.data['label']:
                # get the lineage of the taxid
                lineage = taxa_data.link_taxid(int(taxid),rank)
                if lineage is None:
                    continue
                else:
                    taxa_at_rank.append(lineage)
            taxa_at_rank = np.array(taxa_at_rank)
            for taxa in set(taxa_at_rank):
                taxa_mask = taxa_at_rank == taxa
                taxa_sequences = self.data['seq'][taxa_mask]
                # divide N amongst the taxa sequences
                num_samples = max(1, int(N/sum(taxa_mask)))
                for seq in taxa_sequences:
                    if len(seq) < fragment_size:
                        continue
                    valid_fragments = sample_valid_fragments(seq,num_samples, fragment_size)
                    fragments.extend(valid_fragments)
                    labels.extend([taxa] * len(valid_fragments))

        # overwrite old data with new fragments
        self.data = np.empty(len(labels), dtype=dtype)
        self.data['label'] = labels
        self.data['seq'] = fragments
        self.data['clen'] = np.zeros(len(labels))
    

