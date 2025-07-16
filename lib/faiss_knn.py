import numpy as np
import faiss
import torch

class FaissKNN:
    def __init__(self, dim, metric= 'cosine', use_gpu=True, gpu_device=0):
        """
        Initialize the FAISS KNN index.
        Args:
            dim (int): Dimension of the embeddings.
            metric (str): 'cosine' or 'euclidean'.
            use_gpu (bool): Whether to use GPU (default: True).
            gpu_device (int): GPU device id (default: 0).
        """
        self.dim = dim
        self.metric = metric
        self.use_gpu = use_gpu
        self.gpu_device = gpu_device
        if metric == 'cosine':
            cpu_index = faiss.IndexFlatIP(dim)
            self._normalize = True
        elif metric == 'euclidean':
            cpu_index = faiss.IndexFlatL2(dim)
            self._normalize = False
        else:
            raise ValueError("Unsupported metric. Use 'cosine' or 'euclidean'.")

        if use_gpu:
            res = faiss.StandardGpuResources()
            self.index = faiss.index_cpu_to_gpu(res, gpu_device, cpu_index)
        else:
            self.index = cpu_index

    def add(self, embeddings):
        """
        Add embeddings to the index.
        Args:
            embeddings (np.ndarray): Array of shape (N, D).
        """
        embeddings = np.asarray(embeddings, dtype=np.float32)
        if self._normalize:
            faiss.normalize_L2(embeddings)
        self.index.add(embeddings)

    def search(self, queries, k=5):
        """
        Search for k nearest neighbors.
        Args:
            queries (np.ndarray): Query vectors of shape (M, D).
            k (int): Number of neighbors to return.
        Returns:
            distances (np.ndarray): Distances or similarities of shape (M, k).
            indices (np.ndarray): Indices of neighbors of shape (M, k).
        """
        queries = np.asarray(queries, dtype=np.float32)
        if self._normalize:
            faiss.normalize_L2(queries)
        distances, indices = self.index.search(queries, k)
        return distances, indices 

    def save(self, path):
        """
        Save the FAISS index to disk.
        Args:
            path (str): Path to save the index file.
        """
        # If index is on GPU, move to CPU before saving
        if hasattr(self.index, 'getDevice'):  # GPU index
            cpu_index = faiss.index_gpu_to_cpu(self.index)
            faiss.write_index(cpu_index, path)
        else:
            faiss.write_index(self.index, path)

    @classmethod
    def load(cls, path, metric='cosine', use_gpu=True, gpu_device=0):
        """
        Load a FAISS index from disk.
        Args:
            path (str): Path to the saved index file.
            metric (str): 'cosine' or 'euclidean'.
            use_gpu (bool): Whether to use GPU (default: True).
            gpu_device (int): GPU device id (default: 0).
        Returns:
            FaissKNN instance with loaded index.
        """
        index = faiss.read_index(path)
        if use_gpu:
            res = faiss.StandardGpuResources()
            index = faiss.index_cpu_to_gpu(res, gpu_device, index)
        dim = index.d
        obj = cls(dim=dim, metric=metric, use_gpu=False)  # will overwrite index
        obj.index = index
        return obj

if __name__ == "__main__":
    # Test FaissKNN with random data
    N, D = 100_000, 128
    np.random.seed(42)
    database = np.random.randn(N, D).astype('float32')
    query = database[0].reshape(1, -1)
    k = 10
    use_gpu = torch.cuda.is_available()

    print("--- Cosine Similarity (GPU) ---")
    knn_cosine = FaissKNN(dim=D, metric='cosine', use_gpu=use_gpu)
    knn_cosine.add(database)
    distances, indices = knn_cosine.search(query, k)
    print("Indices:", indices)
    print("Cosine similarities:", distances)

    print("\n--- Euclidean Distance (GPU) ---")
    knn_euclidean = FaissKNN(dim=D, metric='euclidean', use_gpu=use_gpu)
    knn_euclidean.add(database)
    distances, indices = knn_euclidean.search(query, k)
    print("Indices:", indices)
    print("Euclidean distances:", distances)