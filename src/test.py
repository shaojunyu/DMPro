import torch
# pearsonr is a function from scipy.stats
from scipy.stats import pearsonr
import numpy as np

a = np.random.rand(128, 1000)
b = np.random.rand(128, 1000)
# b = np.random.rand(100)
# b = a + np.random.rand(100) * 0.1

# calculate the pearson correlation between a and b, column-wise
pearsonr_list = [pearsonr(a[:, i], b[:, i])[0] for i in range(a.shape[1])]
print(np.mean(pearsonr_list))
print(pearsonr_list[10])
# print(pearsonr(a, b)[0])

# implementaion of pearsonr using pytorch for 2D tensor
def pearsonr_torch(x, y):
    mean_x = torch.mean(x, dim=0)
    mean_y = torch.mean(y, dim=0)
    xm = x - mean_x
    ym = y - mean_y
    
    # dot product of xm and ym in each column
    r_num = xm.t().mm(ym).diagonal()
    r_den = (xm ** 2).sum(dim=0).sqrt() * (ym ** 2).sum(dim=0).sqrt()
    r_val = r_num / r_den
    return r_val


def pearson_loss(y_true, y_pred, eps=1e-8):
    mean_true = torch.mean(y_true)
    mean_pred = torch.mean(y_pred)
    xm = y_true - mean_true
    ym = y_pred - mean_pred
    r_num = xm.dot(ym)
    r_den = (xm ** 2).sum().sqrt() * (ym ** 2).sum().sqrt()
    r_val = r_num / (r_den + eps)
    return (1 - r_val) / 2

a_torch = torch.tensor(a)
b_torch = torch.tensor(b)
res = pearsonr_torch(a_torch, b_torch)
print(res[10])
# print(pearson_loss(a_torch, b_torch))