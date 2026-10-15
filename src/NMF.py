import torch
from torchnmf.nmf import NMF
import numpy as np
import pandas as pd
import argparse
import glob

parser = argparse.ArgumentParser(description='NMF')
parser.add_argument('--rank', type=int, default=400, help='rank')
parser.add_argument('--chr', type=str, help="chromosome", required=True)
parser.add_argument('--encoding_path', type=str, help="path to read the encoding", required=True)
parser.add_argument('--out_path', type=str, help='output file path', required=True)
parser.add_argument('--data_path', type=str, help='data path', required=True)
parser.add_argument('--seed', type=int, default=None,
                    help="random seed for NMF initialisation")
parser.add_argument('--source_platform', type=str, default='HM450',
                    help="platform imputed from, e.g. HM450 or EPIC_V1")
parser.add_argument('--target_platform', type=str, default='EPIC',
                    help="platform imputed to, e.g. EPIC or EPIC_V2")

args = parser.parse_args()

if args.seed is not None:
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

# chr = "chr22"
# encoding_450k = pd.read_csv("../tmp/encoded/HM450_{}_4000.csv".format(chr), index_col=0)
# encoding_EPIC = pd.read_csv("../tmp/encoded/EPIC_{}_8000.csv".format(chr), index_col=0)
# Load data

def _find_encoding(prefix):
    # match the encode output regardless of the bin_size in the filename
    hits = sorted(glob.glob(
        args.encoding_path + "/{}_chr{}_*.csv".format(prefix, args.chr)))
    if not hits:
        raise FileNotFoundError(
            "No encoding for {} chr{} in {}".format(
                prefix, args.chr, args.encoding_path))
    return hits[0]

encoding_src = pd.read_csv(_find_encoding(args.source_platform), index_col=0)
encoding_tgt = pd.read_csv(_find_encoding(args.target_platform), index_col=0)
train_samples = pd.read_csv(
    args.data_path + "/train_samples.txt", header=None
)[0].values
val_samples = pd.read_csv(
    args.data_path + "/val_samples.txt", header=None
)[0].values
test_samples = pd.read_csv(
    args.data_path + "/test_samples.txt", header=None
)[0].values

encoding_src_test = encoding_src.loc[test_samples]
encoding_tgt_test = encoding_tgt.loc[test_samples]

encoding_src_train = encoding_src.drop(test_samples)
encoding_tgt_train = encoding_tgt.drop(test_samples)

# NMF
# concat the source and target encodings of the non-test samples
encoding_train = pd.concat(
    [encoding_src_train, encoding_tgt_train], axis=1).to_numpy()
encoding_train = torch.from_numpy(encoding_train).float()

### first round of NMF
rank = args.rank
model = NMF(encoding_train.shape, rank=rank)
model.fit(encoding_train, tol=1e-10, max_iter=5000, verbose=True)

W = model.W.detach().numpy()
H = model.H.detach().numpy()

# compute the reconstruction error
encoding_train_reconstruct = np.matmul(H, W.T)
encoding_train_reconstruct.shape

np.abs((encoding_train.numpy() - encoding_train_reconstruct)).mean()
# split W at the source column count (the two blocks need not be equal)
n_src = encoding_src.shape[1]
W1 = W[0:n_src, :]
W2 = W[n_src:, :]

# second round of NMF
# source -> target
# keep W1 fixed; torchnmf ignores W= when the data shape is also passed
model_src2tgt = NMF(W=torch.from_numpy(W1).float(),
                    H=(encoding_src_test.shape[0], rank),
                    trainable_W=False)

model_src2tgt.fit(torch.from_numpy(encoding_src_test.to_numpy()),
                  tol=1e-10, max_iter=5000, verbose=True)

# reconstruct the target encoding from the source encoding
encoding_tgt_reconstruct = np.matmul(model_src2tgt.H.detach().numpy(), W2.T)

# save the results
df = pd.DataFrame(encoding_tgt_reconstruct,
                  index=encoding_tgt_test.index, columns=encoding_tgt_test.columns)

print(np.sqrt(((encoding_tgt_test.to_numpy() - encoding_tgt_reconstruct)**2).mean()))

# print(df)
df.to_csv(args.out_path + "/chr{}_reconstruct.csv".format(args.chr))