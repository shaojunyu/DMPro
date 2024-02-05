import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim

data = pd.read_csv('tmp/NMF/EPIC_chr21.csv', index_col=0)

# drop the rows with missing values
data.dropna(inplace=True)
# column-wise normalization, z-score
data = (data - data.mean()) / data.std()

data = data.to_numpy().T
print(data.shape)

# Model Parameters
input_dim = data.shape[1]  # Number of probes
d_model = 512  # Size of the transformer's internal representation
nhead = 8  # Number of heads in multi-head attention
num_encoder_layers = 3
num_decoder_layers = 3
dim_feedforward = 2048

class Model(nn.Module):
    def __init__(self):
        super(Model, self).__init__()
        # encoder layers
        self.encoder_layer = nn.TransformerEncoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=dim_feedforward)
        self.encoder = nn.TransformerEncoder(self.encoder_layer, num_encoder_layers)
        
        # decoder layers
        self.decoder_layer = nn.TransformerDecoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=dim_feedforward)
        self.decoder = nn.TransformerDecoder(self.decoder_layer, num_decoder_layers)
        
        # Input and output linear transformation MLP
        # self.input_linear = nn.Sequential(
        #     nn.Linear(input_dim, 1024),
        #     nn.ReLU(),
        #     nn.Linear(1024, d_model)
        # )
        # self.output_linear = nn.Sequential(
        #     nn.Linear(d_model, 1024),
        #     nn.ReLU(),
        #     nn.Linear(1024, input_dim)
        # )
        self.input_linear = nn.Linear(input_dim, d_model)
        self.output_linear = nn.Linear(d_model, input_dim)
    
    def forward(self, src, tgt):
        src = self.input_linear(src)
        tgt = self.input_linear(tgt)
        memory = self.encoder(src)
        out = self.decoder(tgt, memory)
        out = self.output_linear(out)
        return out

model = Model().cuda()
criterion = nn.L1Loss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

data = torch.tensor(data, dtype=torch.float32) # n x d
data = data.unsqueeze(1) # n x 1 x d
data_loader = torch.utils.data.DataLoader(data, batch_size=32, shuffle=True)
print(data.shape)
# model(data, data)

for epoch in range(1000):
    for batch in data_loader:
        batch = batch.cuda()
        target = batch.clone()
        optimizer.zero_grad()
        output = model(batch, target)
        loss = criterion(output, target)
        loss.backward()
        optimizer.step()
    print('Epoch: {}, Loss: {}'.format(epoch, loss.item()))