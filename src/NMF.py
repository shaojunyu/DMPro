import torch
from torchnmf.nmf import NMF
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='NMF')
parser.add_argument('--rank', type=int, default=400, help='rank')
parser.add_argument('--chr', type=str, help="chromosome", required=True)
parser.add_argument('--encoding_path', type=str, help="path to read the encoding", required=True)
parser.add_argument('--out_path', type=str, help='output file path', required=True)
parser.add_argument('--data_path', type=str, help='data path', required=True)

args = parser.parse_args()

# chr = "chr22"
# encoding_450k = pd.read_csv("../tmp/encoded/HM450_{}_4000.csv".format(chr), index_col=0)
# encoding_EPIC = pd.read_csv("../tmp/encoded/EPIC_{}_8000.csv".format(chr), index_col=0)
# Load data

encoding_450K = pd.read_csv(
    args.encoding_path + "/HM450_chr{}_4000.csv".format(args.chr), index_col=0 
)
encoding_EPIC = pd.read_csv(
    args.encoding_path + "/EPIC_chr{}_8000.csv".format(args.chr), index_col=0
)
train_samples = pd.read_csv(
    args.data_path + "/train_samples.txt", header=None
)[0].values
val_samples = pd.read_csv(
    args.data_path + "/val_samples.txt", header=None
)[0].values
test_samples = pd.read_csv(
    args.data_path + "/test_samples.txt", header=None
)[0].values

encoding_450K_test = encoding_450K.loc[test_samples]
encoding_EPIC_test = encoding_EPIC.loc[test_samples]

encoding_450K_train = encoding_450K.drop(test_samples)
encoding_EPIC_train = encoding_EPIC.drop(test_samples)

# NMF
# concat encoding_450k_train and encoding_EPIC_train
encoding_train = pd.concat(
    [encoding_450K_train, encoding_EPIC_train], axis=1).to_numpy()
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
W1 = W[0:W.shape[0]//2, :]
W2 = W[W.shape[0]//2:, :]

# second round of NMF
# 450k -> EPIC
model_450k2EPIC = NMF(encoding_450K_test.shape, W=W1,
                      trainable_W=False, rank=rank)

model_450k2EPIC.fit(torch.from_numpy(
    encoding_450K_test.to_numpy()), tol=1e-10, max_iter=5000, verbose=True)

# reconstruct EPIC from 450k
encoding_EPIC_reconstruct = np.matmul(model_450k2EPIC.H.detach().numpy(), W2.T)

# save the results
df = pd.DataFrame(encoding_EPIC_reconstruct,
                  index=encoding_EPIC_test.index, columns=encoding_EPIC_test.columns)

print(np.sqrt(((encoding_EPIC_test.to_numpy() - encoding_EPIC_reconstruct)**2).mean()))

# print(df)
df.to_csv(args.out_path + "/chr{}_reconstruct.csv".format(args.chr))