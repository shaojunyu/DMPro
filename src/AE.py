# train an autoencoder for 450K data
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

chr = 15
bin_size = 8000
seed, epochs, batch_size, learning_rate = 0, 500, 128, 0.001
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# Load the data, rows are features, columns are samples
data = pd.read_csv("./tmp/processed/GSE55763_chr{}.csv".format(chr), index_col="ID")
probe_names = data.index.values
sample_names = data.columns.values

# transpose the data so that rows are samples, columns are features
data = data.transpose()

# Handle NA values: Here, I'll replace NA with the mean of the column, but you can choose any other strategy.
data.fillna(data.mean(), inplace=True)

# split the data into training, validation, and testing sets
train_data = data.sample(frac=0.8, random_state=seed)
val_data = data.drop(train_data.index).sample(frac=0.5, random_state=seed)
test_data = data.drop(train_data.index).drop(val_data.index)

# generate probe bins, each bin contains probes of bin_size, bins overlap by bin_size/2 probes
probe_bins = []
for i in range(0, len(probe_names), bin_size // 10):
    if i + bin_size > len(probe_names):
        break
    probe_bins.append(probe_names[i:i + bin_size])

# randomly select some probes, add bins from the selected probes to the training set
# for i in range(0, int(len(probe_names) * 0.5)):
#     probe = np.random.randint(0, len(probe_names))
#     if probe + bin_size > len(probe_names):
#         break
#     probe_bins.append(probe_names[probe:probe + bin_size])

# for each bin, generate a training set
train_sets = []
for bin in probe_bins:
    train_sets.append(train_data[bin].to_numpy())

# test_sets = []
# for bin in probe_bins:
#     test_sets.append(test_data[bin].to_numpy())
val_sets = []
for bin in probe_bins:
    val_sets.append(val_data[bin].to_numpy())
    
# concatenate the training sets, n_samples x bin_size, each row is a training sample
train_data = np.concatenate(train_sets, axis=0)
val_data = np.concatenate(val_sets, axis=0)

print("train_data.shape: ", train_data.shape)

# convert the data to torch tensors, and add a channel dimension, n_samples x 1 x bin_size
train_data = torch.tensor(train_data, dtype=torch.float32).unsqueeze(1).to(device)
val_data = torch.tensor(val_data, dtype=torch.float32).unsqueeze(1).to(device)

# define the dataloader
train_loader = torch.utils.data.DataLoader(train_data, batch_size=batch_size, shuffle=True)
val_loader = torch.utils.data.DataLoader(val_data, batch_size=batch_size, shuffle=False)

class MaskedAutoencoder1D(nn.Module):
    def __init__(self, bin_size):
        super(MaskedAutoencoder1D, self).__init__()
        
        # Compute the spatial size after the convolutional layers
        self.encoded_size = bin_size
        for _ in range(3):  # Three convolutional layers
            self.encoded_size = (self.encoded_size + 2*1 - 3) // 2 + 1
            
        self.encoder = nn.Sequential(
            nn.Conv1d(1, 32, kernel_size=3, stride=2, padding=1),
            nn.ReLU(),
            nn.Conv1d(32, 64, kernel_size=3, stride=2, padding=1),
            nn.ReLU(),
            nn.Conv1d(64, 128, kernel_size=3, stride=2, padding=1),
            nn.Flatten(),
            nn.Linear(128*self.encoded_size, 512),
            nn.ReLU(),
            nn.Linear(512, 128),
            nn.ReLU()
        )
        
        # Decoder
        self.decoder = nn.Sequential(
            nn.Linear(128, 512),
            nn.ReLU(),
            nn.Linear(512, 128*self.encoded_size),
            nn.ReLU(),
            nn.Unflatten(1, (128, self.encoded_size)),
            nn.ConvTranspose1d(128, 64, kernel_size=3, stride=2, padding=1, output_padding=1),
            nn.ReLU(),
            nn.ConvTranspose1d(64, 32, kernel_size=3, stride=2, padding=1, output_padding=1),
            nn.ReLU(),
        )
        
            
        # Encoder
        self.enc1 = nn.Conv1d(1, 32, kernel_size=3, stride=2, padding=1)
        self.enc2 = nn.Conv1d(32, 64, kernel_size=3, stride=2, padding=1)
        self.enc3 = nn.Conv1d(64, 128, kernel_size=3, stride=2, padding=1)
        
        
        # Encoder MLP layers
        self.fc_enc1 = nn.Linear(128*self.encoded_size, 512)
        self.fc_enc2 = nn.Linear(512, 128)
        
        # Decoder
        self.fc_dec1 = nn.Linear(128, 512)
        self.fc_dec2 = nn.Linear(512, 128*self.encoded_size)
        
        self.dec1 = nn.ConvTranspose1d(128, 64, kernel_size=3, stride=2, padding=1, output_padding=1)
        self.dec2 = nn.ConvTranspose1d(64, 32, kernel_size=3, stride=2, padding=1, output_padding=1)
        self.dec3 = nn.ConvTranspose1d(32, 1, kernel_size=3, stride=2, padding=1, output_padding=1)
        
    def forward(self, x, mask_prob=0.2):
        # Apply mask
        mask = (torch.rand_like(x) < mask_prob).float()
        x_masked = x * (1 - mask)
        
        # Encoding
        x = F.relu(self.enc1(x_masked))
        x = F.relu(self.enc2(x))
        x = F.relu(self.enc3(x))
        x = x.view(x.size(0), -1)
        
        # Encoder MLP
        x = F.relu(self.fc_enc1(x))
        encoded = F.relu(self.fc_enc2(x))
        
        # Decoding MLP
        x = F.relu(self.fc_dec1(encoded))
        x = F.relu(self.fc_dec2(x))
        
        x = x.view(x.size(0), 128, self.encoded_size)
        
        # Decoding Conv
        x = F.relu(self.dec1(x))
        x = F.relu(self.dec2(x))
        x = self.dec3(x)
        
        return x


    
# Model, Loss, and Optimizer
model = MaskedAutoencoder1D(bin_size=bin_size).to(device)
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=learning_rate)
scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.1, patience=10, verbose=True)

# Training loop
for epoch in range(epochs):
    model.train()
    train_loss = []
    for data in train_loader:
        # data = data.to(device)
        # Forward pass
        outputs = model(data, mask_prob=0.2)
        loss = criterion(outputs, data)
        train_loss.append(loss.item())
        
        # Backward and optimize
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
    # print('Epoch [{}/{}], Loss: {:.4f}'.format(epoch + 1, epochs, np.mean(train_loss)))
    
    # validation
    val_loss = []
    model.eval()
    with torch.no_grad():
        for data in val_loader:
            # data = data.to(device)
            outputs = model(data)
            loss = criterion(outputs, data)
            val_loss.append(loss.item())
    scheduler.step(np.mean(val_loss))
    
    # print('Epoch [{}/{}], Val Loss: {:.4f}'.format(epoch + 1, epochs, np.mean(val_loss)))
    print("Chr{} Epoch [{}/{}], Train Loss: {:.5f}, Val Loss: {:.5f}".format(chr, epoch + 1, epochs, np.mean(train_loss), np.mean(val_loss)))