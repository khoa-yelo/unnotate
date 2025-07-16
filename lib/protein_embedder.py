"""
Protein Embedding Module using ESM-C
"""

import torch
import pandas as pd
import time
from typing import Dict, Optional, Union, Tuple
import numpy as np
from tqdm import tqdm
def read_fasta_file(fasta_path: str) -> Tuple[list, list]:
    """
    Read protein sequences and IDs from a FASTA file using Biopython.
    
    Args:
        fasta_path (str): Path to the FASTA file.
        
    Returns:
        tuple: (sequences, ids) where sequences is a list of protein sequences and ids is a list of protein IDs.
    """
    try:
        from Bio import SeqIO
    except ImportError as e:
        raise ImportError(f"Biopython not found. Please install it: pip install biopython. Error: {e}")
    
    sequences = []
    ids = []
    
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences.append(str(record.seq))
        ids.append(record.id)
    
    return sequences, ids

class ProteinEmbedder:
    """
    A class for embedding proteins using the ESM-C model.
    """
    
    def __init__(self, model_name: str = "esmc_600m", device: Optional[str] = None):
        """
        Initialize the ProteinEmbedder.
        
        Args:
            model_name (str): Name of the ESM-C model to use.
            device (str, optional): Device to use ('cuda', 'cpu', or None for auto-detect).
        """
        try:
            from esm.models.esmc import ESMC
            from esm.sdk.api import ESMProtein, LogitsConfig
            self.ESMC = ESMC
            self.ESMProtein = ESMProtein
            self.LogitsConfig = LogitsConfig
        except ImportError as e:
            raise ImportError(f"ESM library not found. Please install it: {e}")
        
        # Set device
        if device is None:
            self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        else:
            self.device = torch.device(device)
        
        # Initialize the model
        self.model_name = model_name
        self.client = ESMC.from_pretrained(model_name).to(self.device)
        
        # Default result for empty/invalid sequences
        self.default_result = {
            "max": np.zeros((1152,), dtype=np.float32),
            "mean": np.zeros((1152,), dtype=np.float32),
            "max_middle_layer_12": np.zeros((1152,), dtype=np.float32),
            "mean_middle_layer_12": np.zeros((1152,), dtype=np.float32)
        }
    
    def embed_proteins(self, protein_sequences: Union[list, np.ndarray], 
                      protein_ids: Optional[Union[list, np.ndarray]] = None,
                      middle_layer: int = 12) -> Dict[str, Dict[str, np.ndarray]]:
        """
        Embed multiple protein sequences.
        
        Args:
            protein_sequences (list or np.ndarray): List of protein sequences.
            protein_ids (list or np.ndarray, optional): List of protein IDs. If None, uses indices.
            middle_layer (int): The middle layer to extract logits from.
            
        Returns:
            dict: Dictionary mapping protein IDs to their embeddings.
        """
        if protein_ids is None:
            protein_ids = [f"protein_{i}" for i in range(len(protein_sequences))]
        
        embeddings = {}
        start_time = time.time()
        
        for i, (protein_seq, protein_id) in tqdm(enumerate(zip(protein_sequences, protein_ids)), total=len(protein_sequences), desc="Embedding proteins"):
            # Check if the protein sequence is empty or NaN
            if not protein_seq or pd.isna(protein_seq):
                embeddings[protein_id] = self.default_result.copy()
                continue
            
            protein = self.ESMProtein(sequence=protein_seq)
            protein_tensor = self.client.encode(protein)
            logits_output = self.client.logits(
                protein_tensor,
                self.LogitsConfig(sequence=True, return_embeddings=True, return_hidden_states=True)
            )
            
            embeddings_tensor = logits_output.embeddings.squeeze()
            logits = logits_output.hidden_states[middle_layer].squeeze()
            max_ = torch.max(embeddings_tensor, dim=0).values.float().detach().cpu().numpy()
            mean = embeddings_tensor.mean(dim=0).float().detach().cpu().numpy()
            max_middle = torch.max(logits, dim=0).values.float().detach().cpu().numpy()
            mean_middle = logits.mean(dim=0).float().detach().cpu().numpy()
            embeddings[protein_id] = {
                "max": max_,
                "mean": mean,
                f"max_middle_layer_{middle_layer}": max_middle,
                f"mean_middle_layer_{middle_layer}": mean_middle
            }
        
        return embeddings

# Convenience function for quick embedding
def embed_proteins(protein_sequences: Union[list, str], model_name: str = "esmc_600m", 
                  device: Optional[str] = None, middle_layer: int = 12,
                  protein_ids: Optional[Union[list, np.ndarray]] = None,
                  fasta_file: Optional[str] = None) -> Dict[str, np.ndarray]:
    """
    Convenience function to embed protein sequences.
    
    Args:
        protein_sequences (list or str): List of protein sequences or single protein sequence.
        model_name (str): Name of the ESM-C model to use.
        device (str, optional): Device to use ('cuda', 'cpu', or None for auto-detect).
        middle_layer (int): The middle layer to extract logits from.
        protein_ids (list or np.ndarray, optional): List of protein IDs. If None, uses indices.
        fasta_file (str, optional): Path to FASTA file. If provided, overrides protein_sequences.
        
    Returns:
        dict: A dictionary with metric names as keys and 2D arrays (N, D) as values.
    """
    embedder = ProteinEmbedder(model_name=model_name, device=device)
    
    # Handle FASTA file input
    if fasta_file is not None:
        protein_sequences, protein_ids = read_fasta_file(fasta_file)
    else:
        # Handle single protein sequence
        if isinstance(protein_sequences, str):
            protein_sequences = [protein_sequences]
        
        # Generate IDs if not provided
        if protein_ids is None:
            protein_ids = [str(i) for i in range(len(protein_sequences))]
    
    # Get embeddings from the embedder
    embeddings_dict = embedder.embed_proteins(protein_sequences, protein_ids, middle_layer)
    
    # Reorganize the output format
    N = len(protein_sequences)
    D = 1152  # Embedding dimension
    
    # Initialize output arrays
    result = {
        'max': np.zeros((N, D), dtype=np.float32),
        'mean': np.zeros((N, D), dtype=np.float32),
        f'max_middle_layer_{middle_layer}': np.zeros((N, D), dtype=np.float32),
        f'mean_middle_layer_{middle_layer}': np.zeros((N, D), dtype=np.float32)
    }
    
    # Fill the arrays
    for i, protein_id in enumerate(protein_ids):
        for metric in result.keys():
            result[metric][i] = embeddings_dict[protein_id][metric]
    # sanity check that vectors are not all zeros
    if np.all(result['max'] == 0):
        raise ValueError("All embeddings are zero. Please check your input sequences.")
    return result 