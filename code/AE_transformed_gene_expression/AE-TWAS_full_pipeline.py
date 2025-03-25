"""
AE-TWAS noise robustness tests

This script takes as input the expression modules generated in the previous WGCNA script,
as well as a scale for the *standard devision* of the nonlinear noise to be added. 
Additionally, the script also takes in a suffix that determines whether the noise will be
added to the real data from the previous step, or if a synthetic dataset will be simulated 
based on the real data first before adding the nonlinear noise

By default the intermediate files (including simulated data) generated will be stored in the current working directory, which can be changed via the base_dir variable. 
The input_data_folder, intermediates folder, and model output checkpoints folder will all 
be subdirectories of the base directory. 

"""

# Step 1: Import necessary packages
import pandas as pd
import numpy as np
from torchvision import transforms
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
import pickle
import scipy.stats as stats
from sklearn import preprocessing
import os,sys
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
# Step 2: Set up environment and arguments
np.random.seed(20250409)

# Command-line arguments for module number and noise scale
module_name = sys.argv[1] # String of the module number (e.g: 34, 72)
noise_scale = sys.argv[2]  # Scale of noise to add (e.g: 1e-1, 1.0)
suffix = sys.argv[3] # Either "ramdomAddNonLinearBatchNoise" or "addNonlinearNoise_orig"
base_dir='./' # Top level directory where the data files are
input_data_folder = f"{base_dir}Input_Data" 
output_dir = f"{base_dir}Intermediate_Data"
checkpoints_dir = f"{base_dir}Model_Checkpoints"

suffixes = ["ramdomAddNonLinearBatchNoise", "addNonlinearNoise_orig"]

device = 'cuda' if torch.cuda.is_available() else 'cpu'


# Make output dir if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Step 3: Define utility functions
# Function to add noise to the input DataFrame
def add_noise(input_df, noise_scale="5e-3"):
    """
    Adds noise to the input DataFrame.

    Parameters:
    - input_df: DataFrame to which noise is added.
    - noise_scale: Scale of the noise to be added.

    Returns:
    - output_df: DataFrame with added noise.
    """
    # Add random noise to each feature
    noise_coeff = np.random.normal(loc=0, scale=float(noise_scale), size=input_df.shape[1]) + np.ones(input_df.shape[1])
    noisy_df = input_df.mul(noise_coeff, axis=1)
    
    # Add batch effect noise
    batch_effect_rows = noisy_df.sample(n=200)
    batch_effect_added = batch_effect_rows + 1
    noisy_df.update(batch_effect_added)
    
    # Add additional Gaussian noise
    noise = np.random.normal(loc=0, scale=1, size=input_df.shape)
    residual_noise = pd.DataFrame(noise, index=input_df.index, columns=input_df.columns)
    output_df = noisy_df + residual_noise

    return output_df

def load_and_prepare_simulated_data(
        module_name, 
        noise_scale="1e-3", 
        input_data_folder="Input_Data",
        output_dir="Intermediate_Data",
        suffix = "addNonlinearNoise_orig"
        ):
    # Step 4: Load and preprocess data
    pheno_train_real = pd.read_csv(os.path.join(input_data_folder, f"module_train_genes{module_name}.csv"), sep=",")
    pheno_test_real = pd.read_csv(os.path.join(input_data_folder, f"module_test_genes{module_name}.csv"), sep=",")

    # Rename columns and set index
    pheno_train_real.rename(columns={"Unnamed: 0": "SampleID"}, inplace=True)
    pheno_test_real.rename(columns={"Unnamed: 0": "SampleID"}, inplace=True)
    pheno_train_real.set_index("SampleID", inplace=True)
    pheno_test_real.set_index("SampleID", inplace=True)

    # Concatenate train and test data for simulation
    pheno_train_test_real = pd.concat([pheno_train_real, pheno_test_real], axis=0)

    # Generate simulated data
    pheno_train_test_sim = pd.DataFrame(index=pheno_train_test_real.index, columns=pheno_train_test_real.columns)
    correlation_raw = np.random.normal(0.8, 0.1, size=len(pheno_train_test_sim.columns))
    correlation = [r - 0.5 if r >= 1 else r for r in correlation_raw]

    # Simulate hub genes and correlated genes
    pheno_train_test_sim.iloc[:, 0] = np.random.normal(10, 1, size=pheno_train_test_sim.shape[0])
    pheno_train_test_sim.iloc[:, 1] = 0.8 * pheno_train_test_sim.iloc[:, 0] + np.sqrt(1 - 0.8**2) * np.random.randn(pheno_train_test_sim.shape[0])
    pheno_train_test_sim.iloc[:, 2] = 0.9 * pheno_train_test_sim.iloc[:, 0] + np.sqrt(1 - 0.9**2) * np.random.randn(pheno_train_test_sim.shape[0])

    # Simulate remaining genes based on a randomly selected hub gene
    for col_index in range(3, len(pheno_train_test_sim.columns)):
        hub_gene_index = np.random.randint(0, 3)
        pheno_train_test_sim.iloc[:, col_index] = correlation[col_index] * pheno_train_test_sim.iloc[:, hub_gene_index] + np.sqrt(1 - correlation[col_index]**2) * np.random.randn(pheno_train_test_sim.shape[0])

    # Add noise to simulated data
    pheno_train_test_sim.to_csv(os.path.join(output_dir, f"module_train_test_genes_raw_{module_name}_{noise_scale}.{suffix}.csv"))
    pheno_train_test_sim_addnoise = add_noise(pheno_train_test_sim, float(noise_scale))
    # Save simulated data with noise
    pheno_train_test_sim_addnoise.to_csv(os.path.join(output_dir, f"module_train_test_genes_{module_name}_{noise_scale}.{suffix}.csv"))
    # Prepare data for training and testing
    pheno_train = pheno_train_test_sim_addnoise.loc[pheno_train_real.index].to_numpy()
    pheno_test = pheno_train_test_sim_addnoise.loc[pheno_test_real.index].to_numpy()

    # Normalize data
    scaler = preprocessing.MinMaxScaler()
    pheno_train_norm = scaler.fit_transform(pheno_train)
    pheno_test_norm = scaler.transform(pheno_test)

    return pheno_train_norm, pheno_test_norm

def load_and_prepare_data(
        module_name, 
        noise_scale="1e-3", 
        input_data_folder="Input_Data",
        output_dir="Intermediate_Data",
        suffix = "addNonlinearNoise_orig"
        ):
    """
    Loads and prepares the data by adding noise.

    Parameters:
    - module_number: Module number for which data is loaded.
    - noise_scale: Scale of the noise to be added.
    - data_dir: Directory containing the data files.

    Returns:
    - train_data: Normalized training data.
    - test_data: Normalized test data.
    """
    
    train_filename = f"module_train_genes{module_name}.csv"
    test_filename = f"module_test_genes{module_name}.csv"
    
    # Load datasets
    train_df = pd.read_csv(os.path.join(input_data_folder, train_filename), sep=",")
    test_df = pd.read_csv(os.path.join(input_data_folder, test_filename), sep=",")
    
    # Rename columns and concatenate datasets
    train_df.rename(columns={"Unnamed: 0": "SampleID"}, inplace=True)
    test_df.rename(columns={"Unnamed: 0": "SampleID"}, inplace=True)
    all_samples_df = pd.concat([train_df, test_df])
    all_samples_df.index = all_samples_df["SampleID"]
    all_samples_df = all_samples_df.drop(["SampleID"], axis=1)
    
    # Add noise to the dataset
    noisy_df = add_noise(all_samples_df, noise_scale)
    
    # Save the noisy dataset
    noisy_df.to_csv(os.path.join(output_dir, f"module_genes{module_name}_{noise_scale}.{suffix}.csv"))
    
    # Prepare train and test datasets
    train_data = noisy_df.loc[train_df["SampleID"]].to_numpy()
    test_data = noisy_df.loc[test_df["SampleID"]].to_numpy()
    
    # Normalize the data
    scaler = preprocessing.MinMaxScaler()
    train_data_normalized = scaler.fit_transform(train_data)
    test_data_normalized = scaler.transform(test_data)
    
    return train_data_normalized, test_data_normalized

if suffix == suffixes[0]:
    pheno_train_norm, pheno_test_norm = load_and_prepare_simulated_data(module_name, noise_scale, input_data_folder, output_dir, suffix)
elif suffix == suffixes[1]:
    pheno_train_norm, pheno_test_norm = load_and_prepare_data(module_name, noise_scale, input_data_folder, output_dir, suffix)


# Step 5: Define the dataset class
class InputDataset(Dataset):
    def __init__(self, np_array):
        self.x = np_array
        self.n_sample = np_array.shape[0]

    def __len__(self):
        return self.n_sample

    def __getitem__(self, index):
        return self.x[index]
    
# Step 6: Define the model
class AutoEncoder(nn.Module):
    def __init__(self, input_features):
        super(AutoEncoder, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_features, int(input_features / 2)),
            nn.Sigmoid(),
            nn.Linear(int(input_features / 2), int(input_features / 4)),
            nn.Sigmoid(),
            nn.Linear(int(input_features / 4), int(input_features / 8)),
            nn.Sigmoid()
        )
        self.decoder = nn.Sequential(
            nn.Linear(int(input_features / 8), int(input_features / 4)),
            nn.Sigmoid(),
            nn.Linear(int(input_features / 4), int(input_features / 2)),
            nn.Sigmoid(),
            nn.Linear(int(input_features / 2), input_features),
            nn.Sigmoid()
        )

    def forward(self, x):
        y = self.encoder(x)
        x = self.decoder(y)
        return x, y
    
# Step 7: Prepare data loaders
batch_size = 256  # Example batch size
pheno_train_set = InputDataset(pheno_train_norm)
pheno_train_loader = DataLoader(pheno_train_set, batch_size=batch_size, shuffle=True, num_workers=4)
pheno_test_set = InputDataset(pheno_test_norm)
pheno_test_loader = DataLoader(pheno_test_set, batch_size=batch_size, shuffle=False, num_workers=4)
# Step 8: Set up the model, loss function, optimizer, and scheduler
input_features = pheno_train_norm.shape[1]
model = AutoEncoder(input_features)
loss_fn = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=0.0005)
scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=500, gamma=0.2)

# Step 9: Training and testing the AE
num_epochs = 10001
print_interval = 100
do_train = True
do_test = True
patience = 10
patience_counter = 0
best_test_r2 = float('-inf')
model_save_path = os.path.join(checkpoints_dir, f"{module_name}_{noise_scale}_best_model.L5.{suffix}.pth")
train_test_log_path = os.path.join(checkpoints_dir, f"{module_name}_{noise_scale}_best_model.L5.{suffix}.log")
model_saved = False
with open(train_test_log_path, "w") as train_test_log:
    model.to(device)
    for epoch in range(num_epochs):
        if do_train:
            model.train()
            total_loss = 0.0
            output_train = np.zeros_like(pheno_train_norm)
            input_train = np.zeros_like(pheno_train_norm)
            
            for batch_idx, batch_data in enumerate(pheno_train_loader):
                optimizer.zero_grad()
                batch_data = batch_data.float().to(device)
                output, _ = model(batch_data)
                loss = loss_fn(output, batch_data)
                loss.backward()
                optimizer.step()
                total_loss += loss.item()
                output_train[batch_idx * batch_size:(batch_idx + 1) * batch_size] = output.cpu().detach().numpy()
                input_train[batch_idx * batch_size:(batch_idx + 1) * batch_size] = batch_data.cpu().detach().numpy()
            
            scheduler.step()
            
            if (epoch) % print_interval == 0:
                train_r2 = r2_score(input_train, output_train)
                train_test_log.write(f"Epoch {epoch + 1}, Training Loss: {total_loss}, Training R^2: {train_r2}\n")
                # print(f"Epoch {epoch + 1}, Training Loss: {total_loss}, Training R^2: {train_r2}")

        if do_test:
            if (epoch) % print_interval == 0:
                model.eval()
                total_test_loss = 0.0
                output_test = np.zeros_like(pheno_test_norm)
                input_test = np.zeros_like(pheno_test_norm)
                
                with torch.no_grad():
                    for batch_idx, batch_data in enumerate(pheno_test_loader):
                        batch_data = batch_data.float().to(device)
                        output, _ = model(batch_data)
                        total_test_loss += loss_fn(output, batch_data).item()
                        output_test[batch_idx * batch_size:(batch_idx + 1) * batch_size] = output.cpu().detach().numpy()
                        input_test[batch_idx * batch_size:(batch_idx + 1) * batch_size] = batch_data.cpu().detach().numpy()
                
                test_r2 = r2_score(input_test, output_test)
                train_test_log.write(f"Epoch {epoch + 1}, Test Loss: {total_test_loss}, Test R^2: {test_r2}\n")
                # print(f"Epoch {epoch + 1}, Test Loss: {total_test_loss}, Test R^2: {test_r2}")
                
                
                if test_r2 - best_test_r2 > 0.01:
                    best_test_r2 = test_r2
                    patience_counter = 0
                    if test_r2 > 0.3:
                        torch.save(model.state_dict(), model_save_path)
                        model_saved = True
                        print(f"New best test R^2 = {best_test_r2:.4f}. Model saved.")
                else:
                    patience_counter += 1
                    train_test_log.write(f"No improvement in test R^2 for {patience_counter} check(s).")
                    if patience_counter >= patience or (test_r2 - best_test_r2 < -0.02):
                        print(f"Early stopping triggered at epoch {epoch + 1}. Best test R^2 = {best_test_r2:.4f}")
                        
                        torch.save(model.state_dict(), model_save_path)
                        model_saved = True
                        print(f"New best test R^2 = {test_r2:.4f}. Model saved.")
                        break

if not model_saved:
    print(f"No model found for module {module_name} with noise scale {noise_scale}")
    sys.exit()

# Load train and test IDs
train_sample_ids = pd.read_csv(os.path.join(input_data_folder, f"module_train_genes{module_name}.csv"))["Unnamed: 0"]
test_sample_ids = pd.read_csv(os.path.join(input_data_folder, f"module_test_genes{module_name}.csv"))["Unnamed: 0"]


# Construct file names based on module name and noise level
if suffix == suffixes[0]:
    noised_data_file_name = f"module_train_test_genes_{module_name}_{noise_scale}.{suffix}.csv"
elif suffix == suffixes[1]:
    noised_data_file_name = f"module_genes{module_name}_{noise_scale}.{suffix}.csv"


# Step 10 Load noise data
noise_data = pd.read_csv(os.path.join(output_dir,noised_data_file_name))
noise_data.set_index("SampleID", inplace=True)

# Step 11: Prepare train and test datasets
train_data = noise_data.loc[train_sample_ids]
test_data = noise_data.loc[test_sample_ids]

# Step 12: Normalize the train and test data
scaler = preprocessing.MinMaxScaler()
train_data_normalized = scaler.fit_transform(train_data.to_numpy())
test_data_normalized = scaler.transform(test_data.to_numpy())

# Step 13: Define DataLoader
train_batch_size = train_data_normalized.shape[0]
test_batch_size = test_data_normalized.shape[0]

class GeneExpressionDataset(Dataset):
    def __init__(self, data_array):
        self.data = data_array
        self.num_samples = data_array.shape[0]
    
    def __len__(self):
        return self.num_samples
    
    def __getitem__(self, index):
        return self.data[index]

# Step 14: Define the Autoencoder model
input_features = train_data_normalized.shape[1]
output_features = input_features
hidden_layer_1 = input_features // 2
hidden_layer_2 = input_features // 4
hidden_layer_3 = input_features // 8

class Autoencoder(nn.Module):
    def __init__(self):
        super(Autoencoder, self).__init__()
        # Encoder layers
        self.encoder = nn.Sequential(
            nn.Linear(input_features, hidden_layer_1),
            nn.Sigmoid(),
            nn.Linear(hidden_layer_1, hidden_layer_2),
            nn.Sigmoid(),
            nn.Linear(hidden_layer_2, hidden_layer_3),
            nn.Sigmoid()
        )
        # Decoder layers
        self.decoder = nn.Sequential(
            nn.Linear(hidden_layer_3, hidden_layer_2),
            nn.Sigmoid(),
            nn.Linear(hidden_layer_2, hidden_layer_1),
            nn.Sigmoid(),
            nn.Linear(hidden_layer_1, output_features),
            nn.Sigmoid()
        )
    
    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded, encoded
    
# Step 15: Load the trained model
model_path = f"{module_name}_{noise_scale}_best_model.L5.{suffix}.pth"
model_path = os.path.join(base_dir,"Model_Checkpoints", model_path)
autoencoder_model = Autoencoder()
autoencoder_model.load_state_dict(torch.load(model_path, map_location=device))
autoencoder_model.to(device)

# Step 15.1: Generate a final phenotype file
all_pheno_data = pd.concat([train_data, test_data], axis=0)
all_pheno_data_np_before_norm = all_pheno_data.to_numpy()
all_pheno_data_np = scaler.fit_transform(all_pheno_data_np_before_norm)

# Step 15.2: Setup a DataLoader for the final phenotype file
final_batch_size = all_pheno_data_np.shape[0]
all_pheno_dataset = GeneExpressionDataset(all_pheno_data_np)
all_pheno_loader = DataLoader(all_pheno_dataset, batch_size=final_batch_size, shuffle=False, num_workers=8)

# Step 16: Apply the model to the dataset
autoencoder_model.eval()
mse_loss = nn.MSELoss()  # Mean Squared Error for regression
decoded_pheno = np.zeros_like(all_pheno_data_np)
encoded_code = np.zeros((all_pheno_data_np.shape[0], hidden_layer_3))
total_loss = 0

for batch_idx, batch_data in enumerate(all_pheno_loader):
    batch_tensor = batch_data.float().to(device)
    with torch.no_grad():
        decoded, encoded = autoencoder_model(batch_tensor)
        loss = mse_loss(decoded, batch_tensor)
        total_loss += loss.item()
        decoded_np = decoded.cpu().detach().numpy()
        start_idx = batch_idx * final_batch_size
        end_idx = start_idx + decoded.shape[0]
        decoded_pheno[start_idx:end_idx] = decoded_np
        encoded_code[start_idx:end_idx] = encoded.cpu().detach().numpy()

# Calculate and print the loss and R^2 score
print(f'Loss: {total_loss:.4f}')
r2_score_value = r2_score(all_pheno_data_np, decoded_pheno)
print(f'The average R^2 between y and y_hat for final phenotypes is: {r2_score_value:.4f}')

# Step 17: Output hidden code and decoded phenotype
decoded_pheno_df = pd.DataFrame(decoded_pheno, columns=all_pheno_data.columns, index=all_pheno_data.index)
encoded_code_df = pd.DataFrame(encoded_code, index=all_pheno_data.index)

decoded_pheno_file = f"{module_name}_{noise_scale}_Div8_ID670_output_gcta.{suffix}.txt"
encoded_code_file = f"{module_name}_{noise_scale}_Div8_ID670_code_gcta.{suffix}.txt"
input_pheno_file = f"{module_name}_{noise_scale}_Div8_ID670_input_gcta.{suffix}.txt"

# Save the results to CSV files
decoded_pheno_df.to_csv(os.path.join(output_dir ,decoded_pheno_file), sep=" ", header=None, index=True)
decoded_pheno_df.to_csv(os.path.join(output_dir ,"head_" + decoded_pheno_file), sep=" ", header=True, index=True)
encoded_code_df.to_csv(os.path.join(output_dir ,encoded_code_file), sep=" ", header=None, index=True)
encoded_code_df.to_csv(os.path.join(output_dir ,"head_" + encoded_code_file), sep=" ", header=True, index=True)
all_pheno_data.to_csv(os.path.join(output_dir ,input_pheno_file), sep=" ", header=None, index=True)
all_pheno_data.to_csv(os.path.join(output_dir ,"head_" + input_pheno_file), sep=" ", header=True, index=True)
print("===saved data files===")

# Load raw data based on the module size
# If we're testing simulated data, then module_train_test_genes in Intermediate_Data
# Otherwise it should be the module_genes{module}.csv
if suffix == suffixes[0]:
    raw_data_filepath = os.path.join(output_dir, f"module_train_test_genes_raw_{module_name}_{noise_scale}.{suffix}.csv")
elif suffix == suffixes[1]:
    raw_data_filepath = os.path.join(input_data_folder,f"module_genes{module_name}.csv")
raw_data = pd.read_csv(raw_data_filepath)

# Load data with added noise and AE transformed data
noise_added_data_filepath = f"{output_dir}/head_{module_name}_{noise_scale}_Div8_ID670_input_gcta.{suffix}.txt"
ae_transformed_data_filepath = f"{output_dir}/head_{module_name}_{noise_scale}_Div8_ID670_output_gcta.{suffix}.txt"

noise_added_data = pd.read_csv(noise_added_data_filepath, sep=" ")
ae_transformed_data = pd.read_csv(ae_transformed_data_filepath, sep=" ")

def plot_correlation_bars(raw_data, noise_added_data, ae_transformed_data, module_size, noise_level, suffix, plotdir="plots"):
    """
    Plots a boxplot comparing Pearson correlation coefficients of raw data with noise-added and AE-transformed data.

    Parameters:
    - raw_data: DataFrame containing the original gene expression data.
    - noise_added_data: DataFrame containing the gene expression data with added noise.
    - ae_transformed_data: DataFrame containing the gene expression data after AE transformation.
    - module_size: The size of the gene module.
    - noise_level: The scale of noise added to the data.
    """
    # Lists to store Pearson correlation coefficients
    noise_added_correlations = []
    ae_transformed_correlations = []

    # Calculate Pearson correlation coefficients for each column
    for column in raw_data.columns[1:]:
        noise_corr, _ = pearsonr(noise_added_data[column], raw_data[column])
        ae_corr, _ = pearsonr(ae_transformed_data[column], raw_data[column])
        
        noise_added_correlations.append(abs(noise_corr))
        ae_transformed_correlations.append(abs(ae_corr))

    # Print median correlation coefficients
    median_noise_corr = np.median(noise_added_correlations)
    median_ae_corr = np.median(ae_transformed_correlations)
    print(f"Median Pearson correlation for Module {module_size} with Noise Level {noise_level}:")
    print(f"  Noise Added: {median_noise_corr}, AE Transformed: {median_ae_corr}")

    # Create a DataFrame for plotting
    correlation_data = pd.DataFrame({
        "Noise_v_ori": noise_added_correlations,
        "AE-Transformed_v_ori": ae_transformed_correlations
    })

    correlation_data.to_csv(os.path.join(plotdir,f"{module_name}_{noise_scale}_{suffix}.corr_data.csv"), index=None)

    # Plot the boxplot
    ax = correlation_data.boxplot(grid=False)
    ax.set_ylabel('Pearson r correlation')
    if suffix == suffixes[0]:
        ax.set_title(f"Module {module_size} with simulated Noise Level {noise_level}")
    elif suffix == suffixes[1]:
        ax.set_title(f"Module {module_size} with real Noise Level {noise_level}")
    plt.savefig(os.path.join(base_dir, plotdir , f"correlation_plot_module_{module_size}_noise_{noise_level}_{suffix}.png"))

plot_correlation_bars(raw_data, noise_added_data, ae_transformed_data, module_name, noise_scale, suffix)