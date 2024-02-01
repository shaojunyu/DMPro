import torch
import torch.nn as nn
import torch.nn.functional as F

class MaskedAutoencoder1D(nn.Module):
    def __init__(self, bin_size, encodings_size=128):
        super(MaskedAutoencoder1D, self).__init__()

        self.kernel_size = 3
        # Compute the spatial size after the convolutional layers
        self.bin_size = bin_size
        self.encoded_size = encodings_size
        
        for _ in range(3):  # Three convolutional layers
            self.bin_size = (self.bin_size + 2 - self.kernel_size) // 2 + 1

        # Encoder, with layer normalization
        self.enc1 = nn.Conv1d(
            1, 32, kernel_size=self.kernel_size, stride=2, padding=1)
        self.enc2 = nn.Conv1d(
            32, 64, kernel_size=self.kernel_size, stride=2, padding=1)
        self.enc3 = nn.Conv1d(
            64, 128, kernel_size=self.kernel_size, stride=2, padding=1)

        # Encoder MLP layers
        self.fc_enc1 = nn.Linear(128*self.bin_size, 2048)
        self.fc_enc2 = nn.Linear(2048, encodings_size)

        # Decoder
        self.fc_dec1 = nn.Linear(encodings_size, 2048)
        self.fc_dec2 = nn.Linear(2048, 128*self.bin_size)

        self.dec1 = nn.ConvTranspose1d(
            128, 64, kernel_size=self.kernel_size, stride=2, padding=1, output_padding=1)
        self.dec2 = nn.ConvTranspose1d(
            64, 32, kernel_size=self.kernel_size, stride=2, padding=1, output_padding=1)
        self.dec3 = nn.ConvTranspose1d(
            32, 1, kernel_size=self.kernel_size, stride=2, padding=1, output_padding=1)
    
    def encode(self, x):
        # Encoding
        x = F.relu(self.enc1(x))
        x = F.relu(self.enc2(x))
        x = F.relu(self.enc3(x))
        x = x.view(x.size(0), -1)

        # Encoder MLP
        x = F.relu(self.fc_enc1(x))
        encoded = self.fc_enc2(x)
        
        # make sure the encoded vector is within the range of [0, 1]
        encoded =  torch.layer_norm(encoded, encoded.shape)
        encoded = torch.sigmoid(encoded)
        return encoded

    def decode(self, encoded):
        # Decoding MLP
        x = F.relu(self.fc_dec1(encoded))
        x = F.relu(self.fc_dec2(x))

        x = x.view(x.size(0), 128, self.bin_size)

        # Decoding Conv
        x = F.relu(self.dec1(x))
        x = F.relu(self.dec2(x))
        x = self.dec3(x)
        return x

    def forward(self, x, mask_prob=0):
        # Apply mask
        # mask = (torch.rand_like(x) < mask_prob).float()
        # x_masked = x * (1 - mask)
        
        encoded = self.encode(x)

        x = self.decode(encoded)


        return x
