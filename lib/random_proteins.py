import numpy as np
import random
import torch
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig
import sys
sys.path.insert(0,"/home/jovyan/workspace/unnotate/lib")
from faiss_knn import FaissKNN
import pandas as pd
import h5py

class ProteinSearcher:
    # Constants
    AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
    FAISS_INDEX_PATH = "/home/jovyan/workspace/unnotate/database/mean_middle_layer_12_cosine.faiss.index"
    UNIPROT_DATA_PATH = "/home/jovyan/workspace/unnotate/database/parsed_uniprot_swiss_data.csv"
    EMBEDDINGS_PATH = "/home/jovyan/workspace/unnotate/database/embeddings.h5"

    def __init__(self, device=None):
        self.device = device or torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.esmc_model = None
        self.knn = None

    @staticmethod
    def load_embeddings(embeddings_db, metric):
        """Load embeddings from HDF5 file and select metric."""
        with h5py.File(embeddings_db, "r") as f:
            embeddings = {
                "mean": f["mean"][...],
                "max": f["max"][...],
                "mean_middle_layer_12": f["mean_middle_layer_12"][...],
                "max_middle_layer_12": f["max_middle_layer_12"][...]
            }
            data_ids = f["indices"][...]
        
        if metric not in embeddings:
            raise ValueError(f"Invalid metric: {metric}")
            
        return embeddings[metric], data_ids

    @staticmethod
    def random_protein_sequence(length):
        """Generate a random protein sequence of given length."""
        return ''.join(random.choices(ProteinSearcher.AMINO_ACIDS, k=length))

    @staticmethod
    def generate_random_proteins(n=1000, min_len=40, max_len=1000):
        """Generate n random protein sequences with lengths between min_len and max_len."""
        return [ProteinSearcher.random_protein_sequence(random.randint(min_len, max_len)) 
                for _ in range(n)]

    def initialize_esmc(self):
        """Initialize the ESMC model."""
        if self.esmc_model is None:
            self.esmc_model = ESMC.from_pretrained("esmc_600m").to(self.device)

    def embed_proteins(self, seqs, middle_layer=12):
        """Embed a batch of protein sequences with ESMC."""
        self.initialize_esmc()
        embeddings = []
        for seq in seqs:
            protein = ESMProtein(sequence=seq)
            protein_tensor = self.esmc_model.encode(protein)
            logits_output = self.esmc_model.logits(
                protein_tensor,
                LogitsConfig(sequence=True, return_embeddings=True, return_hidden_states=True)
            )
            mean_middle = logits_output.hidden_states[middle_layer].squeeze().mean(dim=0).float().detach().cpu().numpy()
            embeddings.append(mean_middle)
        return np.stack(embeddings)

    def initialize_faiss(self):
        """Initialize the FAISS index."""
        if self.knn is None:
            self.knn = FaissKNN.load(
                self.FAISS_INDEX_PATH, 
                metric='cosine', 
                use_gpu=torch.cuda.is_available()
            )

    def search_proteins(self, proteins, k=40):
        """Search protein embeddings against FAISS index."""
        print("Embedding proteins with ESMC...")
        embeddings = self.embed_proteins(proteins)

        print("Searching against FAISS index...")
        self.initialize_faiss()
        distances, indices = self.knn.search(embeddings, k=k)

        return distances, indices

    def process_results(self, proteins, distances, indices):
        """Process search results and save to file."""
        reference_embeddings, data_ids = self.load_embeddings(self.EMBEDDINGS_PATH, metric='mean_middle_layer_12')
        accession_matrix = data_ids[indices].astype(str)

        df_uniprot = pd.read_csv(self.UNIPROT_DATA_PATH)
        df_accession = df_uniprot[df_uniprot['accession'].isin(accession_matrix.flatten())]
        
        query_lengths = [len(p) for p in proteins]
        accession_length_map = dict(zip(df_accession['accession'].values, df_accession['sequence_length'].values))
        accession_lengths = np.array([accession_length_map[acc] for acc in accession_matrix.flatten()])
        accession_lengths = accession_lengths.reshape(accession_matrix.shape)

        np.savez(
            "calibrate.npz", 
            distances=distances, 
            query_lengths=query_lengths, 
            accession_matrix=accession_matrix, 
            accession_lengths=accession_lengths
        )

def main():
    searcher = ProteinSearcher()
    
    print("Generating random protein sequences...")
    proteins = searcher.generate_random_proteins(n=10000, min_len=30, max_len=1200)
    
    distances, indices = searcher.search_proteins(proteins)
    
    print("Processing and saving results...")
    searcher.process_results(proteins, distances, indices)

if __name__ == "__main__":
    main()