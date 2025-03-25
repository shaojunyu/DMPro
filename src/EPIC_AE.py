# train an autoencoder for 450K data
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import argparse
from model import MaskedAutoencoder1D

parser = argparse.ArgumentParser()
parser.add_argument("--chr", type=str, default="20", help="chromosome number")
parser.add_argument("--bin_size", type=int, default=4000, help="bin size")
parser.add_argument("--seed", type=int, default=0, help="random seed")
parser.add_argument("--epochs", type=int, default=300, help="number of epochs")
parser.add_argument("--batch_size", type=int, default=128, help="batch size")
parser.add_argument("--learning_rate", type=float,
                    default=0.001, help="learning rate")
parser.add_argument("--mode", type=str, default="train",
                    help="train, encode, or test")
parser.add_argument("--model_path", type=str, default="", required=True,
                    help="model path to load or save")
parser.add_argument("--data_prefix", type=str, default="tmp/processed/",
                    required=False, help="data path prefix")
parser.add_argument("--platform", type=str, default="HM450",
                    help="HM450 or EPIC, Illumina platform")
parser.add_argument("--step_size", type=int, default=2000, help="step size")
parser.add_argument("--bin_size_file", type=str, 
                    # default="src/training_bin_size.csv",
                    help="file to store the bin sizes")
parser.add_argument("--output_path", type=str, default="tmp/encoded/",
                    help="output path for encoded data")
parser.add_argument("--encoding_file", type=str, default="",
                    help="file to store the encoded data")
parser.add_argument("--decode_file", type=str, default="",
                    help="file to store the decoded data")
parser.add_argument("--pearson_loss", type=float, default=0,
                    help="ratio of pearson loss")
                    
args = parser.parse_args()

# read the bin size file
if args.bin_size_file:
    bin_size_df = pd.read_csv(args.bin_size_file)
    bin_size_df = bin_size_df[(bin_size_df["platform"] == args.platform) & (bin_size_df["chr"] == int(args.chr))]
    args.bin_size = bin_size_df.bin_size.to_list()[0]
    args.step_size = bin_size_df.step_size.to_list()[0]
    print("bin_size: ", args.bin_size, "step_size: ", args.step_size)

bin_size = args.bin_size
seed, epochs, batch_size, learning_rate = 0, args.epochs, args.batch_size, args.learning_rate
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
# device = torch.device("cpu")

# tmp/processed/EPIC_chr10.csv
# Load the data, rows are features, columns are samples
data_path = args.data_prefix + "/{}_chr{}.csv".format(args.platform, args.chr)
print("Data path: ", data_path)
data = pd.read_csv(data_path, index_col="ID")
data = data.drop(["chr", "MAPINFO"], axis=1)

if args.platform == "EPIC":
    # read the corresponding 450K data, use only the EPIC-specific probes
    # data_450k = pd.read_csv(data_path.replace("EPIC", "HM450"), index_col="ID")
    # keep = data.index.difference(data_450k.index)
    # print("Original data shape: ", data.shape)
    # data = data.loc[keep]
    # print("New data shape: ", data.shape)
    pass

probe_names = data.index.values
print("probe_names: ", len(probe_names))
print("unique probes: ", len(np.unique(probe_names)))
sample_names = data.columns.values

# transpose the data so that rows are samples, columns are features
data = data.transpose()


if data.isna().sum().sum() > 0:
    print("NA values found, fill with mean")
    data.fillna(data.mean(), inplace=True)
# Handle NA values: Here, I'll replace NA with the mean of the column, but you can choose any other strategy.
# data.fillna(data.mean(), inplace=True)
# print("device: ", device)
# print("device: ", data)

# split the data into training, validation, and testing sets
train_samples = pd.read_csv(args.data_prefix + 
    "/train_samples.txt", header=None)[0].values
val_samples = pd.read_csv(args.data_prefix +
    "/val_samples.txt", header=None)[0].values
test_samples = pd.read_csv(args.data_prefix +
    "/test_samples.txt", header=None)[0].values
train_data = data.loc[train_samples]
val_data = data.loc[val_samples]
test_data = data.loc[test_samples]

# train_data = data.sample(frac=0.8, random_state=seed)
# val_data = data.drop(train_data.index).sample(frac=0.5, random_state=seed)
# test_data = data.drop(train_data.index).drop(val_data.index)

if bin_size > len(probe_names):
    bin_size = len(probe_names)
    train_data = train_data.to_numpy()
    val_data = val_data.to_numpy()
    test_data = test_data.to_numpy()

else:
    # generate probe bins, each bin contains probes of bin_size, bins overlap by bin_size/2 probes
    probe_bins = []
    for i in range(0, len(probe_names), args.step_size):
        if i + bin_size > len(probe_names):
            probe_bins.append(probe_names[-bin_size:])
            break
        probe_bins.append(probe_names[i:i + bin_size])

    # for each bin, generate a training set
    train_sets = []
    for bin in probe_bins:
        train_sets.append(train_data[bin].to_numpy())

    val_sets = []
    for bin in probe_bins:
        val_sets.append(val_data[bin].to_numpy())

    # concatenate the training sets, n_samples x bin_size, each row is a training sample
    train_data = np.concatenate(train_sets, axis=0)
    val_data = np.concatenate(val_sets, axis=0)

    # print("train_data.shape: ", train_data.shape)

# convert the data to torch tensors, and add a channel dimension, n_samples x 1 x bin_size
train_data = torch.tensor(
    train_data, dtype=torch.float32).unsqueeze(1).to(device)
val_data = torch.tensor(val_data, dtype=torch.float32).unsqueeze(1).to(device)

print("train_data.shape: ", train_data.shape)

# define the dataloader
train_loader = torch.utils.data.DataLoader(
    train_data, batch_size=batch_size, shuffle=True)
val_loader = torch.utils.data.DataLoader(
    val_data, batch_size=batch_size, shuffle=False)


# Model, Loss, and Optimizer
model = MaskedAutoencoder1D(bin_size=bin_size).to(device)
# criterion = nn.MSELoss()

# use MAE loss to avoid outliers
criterion = nn.L1Loss()

# Loss Function with Pearson’s Correlation Coefficient
def pearson_loss(y_true, y_pred, eps=1e-8):
    '''
    y_true: true values, shape: (n_samples, 1, bin_size)
    y_pred: predicted values, shape: (n_samples, 1, bin_size)
    eps: epsilon, to avoid division by zero
    calculate the Pearson correlation coefficient loss between y_true and y_pred for each probe
    '''
    # reduce the dimension of y_true and y_pred, shape: (n_samples, bin_size)
    y_true = y_true.squeeze(1)
    y_pred = y_pred.squeeze(1)
    mean_x = torch.mean(y_true, dim=0)
    mean_y = torch.mean(y_pred, dim=0)
    xm = y_true - mean_x
    ym = y_pred - mean_y
    # dot product of xm and ym in each column
    r_num = xm.t().mm(ym).diagonal()
    r_den = (xm ** 2).sum(dim=0).sqrt() * (ym ** 2).sum(dim=0).sqrt() + eps
    r_val = r_num / r_den
    return (1 - r_val) / 2

optimizer = optim.Adam(model.parameters(), lr=learning_rate)
scheduler = optim.lr_scheduler.ReduceLROnPlateau(
    optimizer, mode='min', factor=0.1, patience=20, verbose=True)

if args.mode == "encode":
    # model.load_state_dict(torch.load(
    #     "./tmp/models/{}_chr{}_{}.pt".format(args.platform, args.chr, bin_size)))
    model.load_state_dict(torch.load(
        args.model_path + "/{}_chr{}_{}_{}.pt".format(args.platform, args.chr, bin_size, args.step_size)))
    model.eval()
    print("Encoding...")
    # data = data.to_numpy()
    # data = torch.tensor(data, dtype=torch.float32).unsqueeze(1).to(device)
    print("data", data.shape)
    # define new probe bins for testing, increase the stride to reduce the number of bins
    probe_bins = []
    for i in range(0, len(probe_names), args.step_size):
        if i + bin_size > len(probe_names):
            probe_bins.append(probe_names[-bin_size:])
            break
        probe_bins.append(probe_names[i:i + bin_size])
    print("probe_bins", len(probe_bins))

    encoding_sets = []
    for bin in probe_bins:
        current_bin = data[bin]

        # current_bin["sample"] = current_bin.index
        sample_names = current_bin.index.values
        probe_names = current_bin.columns.values

        # encoding
        with torch.no_grad():
            input = torch.tensor(current_bin.to_numpy(),
                                 dtype=torch.float32).unsqueeze(1).to(device)
            encoding = model.encode(input).cpu().numpy()  # n_samples x 128
            encoding_sets.append(encoding)

    encoding = np.concatenate(encoding_sets, axis=1)
    encoding = pd.DataFrame(encoding, index=sample_names)
    print("encoding", encoding.shape)
    # encoding.to_csv(
    #     "./tmp/encoded/{}_chr{}_{}.csv".format(args.platform, args.chr, bin_size))
    # encoding.to_csv()
    encoding.to_csv(
        args.output_path + "/{}_chr{}_{}.csv".format(args.platform, args.chr, bin_size))
    exit()

if args.mode == "decode":
    model.load_state_dict(torch.load(
        args.model_path + "/{}_chr{}_{}_{}.pt".format(args.platform, args.chr, bin_size, args.step_size)))
    model.eval()
    print("Decoding...")

    # encoding = pd.read_csv(
    #     "./tmp/encoded/EPIC_chr{}_{}_reconstruct.csv".format(args.chr, args.bin_size), index_col=0)
    encoding = pd.read_csv(
        args.encoding_file, index_col=0)
        
    sample_names = encoding.index.values
    encoding = encoding.to_numpy()
    encoding = torch.tensor(encoding, dtype=torch.float32).unsqueeze(
        1).to(device)  # n_samples x 1 x 640
    print("encoding", encoding.shape)
    # define new probe bins for testing, increase the stride to reduce the number of bins
    probe_bins = []
    for i in range(0, len(probe_names), args.step_size):
        if i + bin_size > len(probe_names):
            probe_bins.append(probe_names[-bin_size:])
            break
        probe_bins.append(probe_names[i:i + bin_size])

    decode_sets = []
    print("probe_bins", len(probe_bins), encoding.shape[2] / 128)
    for i, bin in enumerate(probe_bins):
        current_bin_probes = bin
        # print(bin)
        decode = model.decode(encoding[:, :, i*128:(i+1)*128])
        decode = decode.squeeze(1)
        decode = pd.DataFrame(decode.detach().cpu().numpy(), index=sample_names,
                              columns=current_bin_probes)
        # print("decode", decode)
        decode["sample"] = decode.index
        decode = pd.melt(decode, id_vars=[
                         "sample"], var_name="probe", value_vars=current_bin_probes)
        decode_sets.append(decode)

    decode = pd.concat(decode_sets, axis=0)
    # group by sample and probe, and take the mean of the predictions
    decode = decode.groupby(["sample", "probe"]).mean()
    decode = decode.reset_index()
    decode = pd.pivot(decode, index="sample", columns="probe", values="value")
    print("decode", decode.shape)
    # transpose the data so that rows are probes, columns are samples
    decode = decode.transpose()
    decode.to_csv(
        args.decode_file)
    # decode.to_csv(
    #     "./tmp/decoded/{}_chr{}_{}.csv".format(args.platform, args.chr, bin_size))
    exit()

if args.mode == "test":
    # define new probe bins for testing, increase the stride to reduce the number of bins
    probe_bins = []
    for i in range(0, len(probe_names), args.step_size):
        if i + bin_size > len(probe_names):
            probe_bins.append(probe_names[-bin_size:])
            break
        probe_bins.append(probe_names[i:i + bin_size])
    print("chr:", args.chr, "probe_bins", len(probe_bins))
    # import glob
    # m = glob.glob(args.model_path + "/{}_chr{}_{}*".format(args.platform,
    #           args.chr, bin_size))
    # exit()
    m = args.model_path + "/{}_chr{}_{}_{}.pt".format(args.platform, args.chr, bin_size, args.step_size)
    model.load_state_dict(torch.load(m))
    # print("Model: ", m[0])
    # load the model
    # model.load_state_dict(torch.load(
    #     args.model_path + "/{}_chr{}_{}_{}.pt".format(args.platform, args.chr, bin_size, args.step_size)))
    # print(model.device)
    # torch.load(m[0], map_location=torch.device('cpu'))
    
    # model.load_state_dict(torch.load(m[0], map_location=device))
    model.eval()

    # print("Testing...")
    # print("Test data", test_data.shape)
    prediction_sets = []
    # merge all probes into one bin
    all_probes = []
    for bin in probe_bins:
        all_probes.extend(bin)
    # print("all_probes", len(all_probes))

    for bin in probe_bins:
        current_bin = test_data[bin]

        # current_bin["sample"] = current_bin.index
        sample_names = current_bin.index.values
        probe_names = current_bin.columns.values
        # pd.melt(current_bin, id_vars=["sample"], var_name="probe", value_vars=probe_names)

        # predict
        with torch.no_grad():
            input = torch.tensor(current_bin.to_numpy(),
                                 dtype=torch.float32).unsqueeze(1).to(device)
            outputs = model(input).cpu().numpy()  # n_samples x 1 x bin_size
            outputs = outputs.squeeze(1)  # n_samples x bin_size
            outputs = pd.DataFrame(
                outputs, index=sample_names, columns=probe_names)
            outputs["sample"] = outputs.index
            outputs = pd.melt(outputs, id_vars=[
                              "sample"], var_name="probe", value_vars=probe_names)
            prediction_sets.append(outputs)

    prediction = pd.concat(prediction_sets, axis=0)
    # group by sample and probe, and take the mean of the predictions
    prediction = prediction.groupby(["sample", "probe"]).mean()
    
    # distinct by sample and probe, keep the first value
    # prediction = prediction.drop_duplicates(subset=["sample", "probe"])
    prediction = prediction.reset_index()
    prediction = pd.pivot(prediction, index="sample",
                          columns="probe", values="value")
    # print("prediction", prediction.shape)

    # rearrange the rows and columns to match the test data
    prediction = prediction.loc[test_data.index.values,
                                test_data.columns.values]
    
    # for each column of prediction and test data, calculate the correlation
    corr = []
    for i in range(len(prediction.columns)):
        corr.append(np.corrcoef(
            prediction.iloc[:, i], test_data.iloc[:, i])[0, 1])
    print("Average correlation: {:.4f}".format(np.mean(corr)))
    
    # for each row of prediction and test data, calculate the correlation and RMSE
    rmse = []
    for i in range(len(prediction)):
        # print("Correlation of row {}: {:.4f}".format(
        #     i, np.corrcoef(prediction.iloc[i], test_data.iloc[i])[0, 1]))
        # print("RMSE of row {}: {:.4f}".format(
        #     i, np.sqrt(np.mean((prediction.iloc[i] - test_data.iloc[i])**2))))
        rmse.append(
            np.sqrt(np.mean((prediction.iloc[i] - test_data.iloc[i])**2)))
    print("Average RMSE: {:.4f}".format(np.mean(rmse)))
    exit()

# Training loop
best_val_loss = np.inf
best_model = None

for epoch in range(epochs):
    model.train()
    train_loss = []
    for data in train_loader:
        # data = data.to(device)
        # Forward pass
        outputs = model(data, mask_prob=0)
        loss = criterion(outputs, data)
        if args.pearson_loss:
            p_loss = pearson_loss(outputs, data)
            loss = loss + p_loss.mean()
        
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
            p_loss = pearson_loss(outputs, data)
            p_loss = p_loss.mean() * args.pearson_loss
            loss = loss + p_loss.mean()
            
            val_loss.append(loss.item())
    scheduler.step(np.mean(val_loss))

    if np.mean(val_loss) < best_val_loss:
        best_val_loss = np.mean(val_loss)
        best_model = model.state_dict()
        # torch.save(model.state_dict(
        # ), args.model_path + "/{}_chr{}_{}_{}.pt".format(args.platform, args.chr, bin_size, args.step_size))

    # compute the RMSE
    # rmse = np.sqrt(np.mean(val_loss**2))

    # print('Epoch [{}/{}], Val Loss: {:.4f}'.format(epoch + 1, epochs, np.mean(val_loss)))
    print("{}, Chr{} Epoch [{}/{}], Train Loss: {:.5f}, Val Loss: {:.5f}".format(
        args.platform, args.chr, epoch + 1, epochs, np.mean(train_loss), np.mean(val_loss)), flush=True)

# save the model
torch.save(best_model, args.model_path +
           "/{}_chr{}_{}_{}.pt".format(args.platform, args.chr, bin_size, args.step_size))
