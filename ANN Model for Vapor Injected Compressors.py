"""
Developed by AMJID KHAN
PhD Student of Dr. Craig Bradshaw
ANN Model for Performance Prediction of Vapor Injected Compressor
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from sklearn.metrics import mean_absolute_percentage_error, mean_absolute_error
from matplotlib.lines import Line2D


"Compressor Power"


plt.figure(figsize=(6,4))
font = {'family': 'Times New Roman', 'color': 'black', 'weight': 'normal', 'size': 14}

# Read data from both sheets
df_training = pd.read_excel('ANN DATA.XLSX', sheet_name=0)
df_testing = pd.read_excel('ANN DATA.XLSX', sheet_name=1)

# Combine data from both sheets
df_combined = pd.concat([df_training, df_testing])

# Split combined data into features and labels
features = df_combined[['P_suc', 'P_dis', 'P_inj', 'F']].to_numpy()

"""
F = compressor Speed 
P_suc = Suction Pressure
P_inj = Injection pressure
P_dis = Discharge Pressure
W_dot_comp = Compressor Power"""


labels = df_combined[['W_dot_comp']].to_numpy()

# Scale features and labels
x_scaler = preprocessing.MinMaxScaler()
y_scaler = preprocessing.MinMaxScaler()

X = x_scaler.fit_transform(features)
y = y_scaler.fit_transform(labels)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.8, random_state=1)

# Train the model
neural_network = MLPRegressor(random_state=1, max_iter=3000, hidden_layer_sizes=(20,), activation='relu', solver='lbfgs')
neural_network.fit(X_train, y_train)

# Predictions
y_pred_test = neural_network.predict(X_test)
y_pred_train = neural_network.predict(X_train)

# Inverse scaling for predictions
y_pred_test = y_scaler.inverse_transform(y_pred_test.reshape(-1, 1))
y_pred_train = y_scaler.inverse_transform(y_pred_train.reshape(-1, 1))
y_test = y_scaler.inverse_transform(y_test)
y_train = y_scaler.inverse_transform(y_train)

# Evaluation
test_error = np.sqrt(mean_squared_error(y_test, y_pred_test))
train_error = np.sqrt(mean_squared_error(y_train, y_pred_train))

print("Test RMSE:", test_error)
print("Train RMSE:", train_error)

# Calculate absolute percentage error
abs_percentage_error_test = np.abs((y_test - y_pred_test) / y_test) * 100

# Calculate mean absolute percentage error
mape_Power = np.mean(abs_percentage_error_test)

print("Mean Absolute Percentage Error (Test):", mape_Power)


# Parity plot for training and test data
plot_limit_min = (df_testing["W_dot_comp"]).min()-200
plot_limit_max = (df_testing["W_dot_comp"]).max()+200
x = np.linspace(plot_limit_min, plot_limit_max, 10)
y = x
y_plus_10 = y * 1.1
y_minus_10 = y * 0.9

plt.scatter(y_train, y_pred_train, alpha=0.9, edgecolor='darkgreen', facecolor='none', marker='^', s=100, linewidths=3, label='Training data')
plt.scatter(y_test, y_pred_test, alpha=0.9, edgecolor='darkblue', facecolor='none', marker='D', s=80, linewidths=2, label='Testing data')

plt.minorticks_on()
plt.tick_params(direction='in', which='minor',bottom=True, top=True, left=True, right=True, length= 3)
plt.tick_params(direction='in', which='major',bottom=True, top=True, left=True, right=True, length= 7)

#plt.colorbar(label="Relative Percent Error",orientation="vertical")
plt.plot(x, y, color='r', alpha=0.4)
plt.plot(x, y_plus_10, linestyle=':', color='k', alpha=0.9, label= " $\pm$ 10 %")
plt.plot(x, y_minus_10, linestyle=':', color='k', alpha=0.9)
plt.xlim(plot_limit_min, plot_limit_max)
plt.ylim(plot_limit_min, plot_limit_max)

plt.xlabel("Actual Compressor Power [W]", fontdict=font)
plt.ylabel("Predicted Compressor Power [W]", fontdict=font)
plt.title("Compressor Power [W]", fontdict=font)

plt.fill_between(x, y, y_plus_10, color='blue', alpha=0.1)

plt.fill_between(x, y, y_minus_10, color='blue', alpha=0.1)

legend_elements = [Line2D([0], [0], marker='D', color='darkblue', linestyle='',
                                          label='Training Data',
                                         markerfacecolor='none', markersize=10),
                   Line2D([0], [0], marker='^', color='darkgreen', linestyle='',
                          label='Testing Data',
                          markerfacecolor='none', markersize=10),
                   Line2D([0], [0], marker='', color='darkgreen', linestyle='',
                          label='MAPE = {:.2f} %'.format(mape_Power),
                          markerfacecolor='none', markersize=10)]
second_legend = plt.legend(handles=legend_elements,loc='upper left', prop={'family': 'Times New Roman', 'size': 12})


plt.tight_layout() 
plt.savefig('Compressor Power', dpi=600)
plt.show() 

#%%


"Evaporator Mass Flow Rate"

plt.figure(figsize=(6,4))
font = {'family': 'Times New Roman', 'color': 'black', 'weight': 'normal', 'size': 14}

# Read data from both sheets
df_training = pd.read_excel('ANN DATA.XLSX', sheet_name=0)
df_testing = pd.read_excel('ANN DATA.XLSX', sheet_name=1)

# Combine data from both sheets
df_combined = pd.concat([df_training, df_testing])

# Split combined data into features and labels
features = df_combined[['P_suc', 'P_dis', 'P_inj', 'F']].to_numpy()

"""
F = compressor Speed 
P_suc = Suction Pressure
P_inj = Injection pressure
P_dis = Discharge Pressure
m_dot_suc = Suction Mass Flow Rate""" 

labels = df_combined[['m_dot_suc']].to_numpy()

# Scale features and labels
x_scaler = preprocessing.MinMaxScaler()
y_scaler = preprocessing.MinMaxScaler()

X = x_scaler.fit_transform(features)
y = y_scaler.fit_transform(labels)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.8, random_state=1)

# Train the model
neural_network = MLPRegressor(random_state=1, max_iter=6000, hidden_layer_sizes=(20,), activation='relu', solver='lbfgs')
neural_network.fit(X_train, y_train)

# Predictions
y_pred_test = neural_network.predict(X_test)
y_pred_train = neural_network.predict(X_train)

# Inverse scaling for predictions
y_pred_test = y_scaler.inverse_transform(y_pred_test.reshape(-1, 1))
y_pred_train = y_scaler.inverse_transform(y_pred_train.reshape(-1, 1))
y_test = y_scaler.inverse_transform(y_test)
y_train = y_scaler.inverse_transform(y_train)

# Evaluation
test_error = np.sqrt(mean_squared_error(y_test, y_pred_test))
train_error = np.sqrt(mean_squared_error(y_train, y_pred_train))

print("Test RMSE:", test_error)
print("Train RMSE:", train_error)
# Calculate absolute percentage error
abs_percentage_error_test = np.abs((y_test - y_pred_test) / y_test) * 100

# Calculate mean absolute percentage error
mape_m_dot_suc = np.mean(abs_percentage_error_test)

print("Mean Absolute Percentage Error (Test):", mape_m_dot_suc)


# Parity plot for training and test data
plot_limit_min = (df_testing["m_dot_suc"]).min()-2
plot_limit_max = (df_testing["m_dot_suc"]).max()+2
x = np.linspace(plot_limit_min, plot_limit_max, 10)
y = x
y_plus_10 = y * 1.1
y_minus_10 = y * 0.9

plt.scatter(y_train, y_pred_train, alpha=0.9, edgecolor='darkgreen', facecolor='none', marker='^', s=100, linewidths=3, label='Training data')
plt.scatter(y_test, y_pred_test, alpha=0.9, edgecolor='darkblue', facecolor='none', marker='D', s=80, linewidths=2, label='Testing data')

plt.minorticks_on()
plt.tick_params(direction='in', which='minor',bottom=True, top=True, left=True, right=True, length= 3)
plt.tick_params(direction='in', which='major',bottom=True, top=True, left=True, right=True, length= 7)

#plt.colorbar(label="Relative Percent Error",orientation="vertical")
plt.plot(x, y, color='r', alpha=0.4)
plt.plot(x, y_plus_10, linestyle=':', color='k', alpha=0.9, label= " $\pm$ 10 %")
plt.plot(x, y_minus_10, linestyle=':', color='k', alpha=0.9)
plt.xlim(plot_limit_min, plot_limit_max)
plt.ylim(plot_limit_min, plot_limit_max)

plt.xlabel("Actual Evaporator Mass Flow Rate [lb/min]", fontdict=font)
plt.ylabel("Predicted Evaporator Mass Flow Rate [lb/min]", fontdict=font)
plt.title("Evaporator Mass Flow Rate [lb/min]", fontdict=font)

plt.fill_between(x, y, y_plus_10, color='blue', alpha=0.1)

plt.fill_between(x, y, y_minus_10, color='blue', alpha=0.1)


legend_elements = [Line2D([0], [0], marker='D', color='darkblue', linestyle='',
                                          label='Training Data',
                                         markerfacecolor='none', markersize=10),
                   Line2D([0], [0], marker='^', color='darkgreen', linestyle='',
                          label='Testing Data',
                          markerfacecolor='none', markersize=10),
                   Line2D([0], [0], marker='', color='darkgreen', linestyle='',
                          label='MAPE = {:.2f} %'.format(mape_Power),
                          markerfacecolor='none', markersize=10)]
second_legend = plt.legend(handles=legend_elements,loc='upper left', prop={'family': 'Times New Roman', 'size': 12})


plt.tight_layout() 
plt.savefig('Suction Mass Flow Rate', dpi=600)
plt.show() 

#%%


"Injection Mass Flow Rate"

plt.figure(figsize=(6,4))
font = {'family': 'Times New Roman', 'color': 'black', 'weight': 'normal', 'size': 14}

# Read data from both sheets
df_training = pd.read_excel('ANN DATA.XLSX', sheet_name=0)
df_testing = pd.read_excel('ANN DATA.XLSX', sheet_name=1)

# Combine data from both sheets
df_combined = pd.concat([df_training, df_testing])

# Split combined data into features and labels
features = df_combined[['P_suc', 'P_dis', 'P_inj', 'F']].to_numpy()

"""
F = compressor Speed 
P_suc = Suction Pressure
P_inj = Injection pressure
P_dis = Discharge Pressure
m_dot_inj = Injection Mass Flow Rate"""

labels = df_combined[['m_dot_inj']].to_numpy()

# Scale features and labels
x_scaler = preprocessing.MinMaxScaler()
y_scaler = preprocessing.MinMaxScaler()

X = x_scaler.fit_transform(features)
y = y_scaler.fit_transform(labels)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.8, random_state=1)

# Train the model
neural_network = MLPRegressor(random_state=1, max_iter=6000, hidden_layer_sizes=(20,), activation='relu', solver='lbfgs')
neural_network.fit(X_train, y_train)

# Predictions
y_pred_test = neural_network.predict(X_test)
y_pred_train = neural_network.predict(X_train)

# Inverse scaling for predictions
y_pred_test = y_scaler.inverse_transform(y_pred_test.reshape(-1, 1))
y_pred_train = y_scaler.inverse_transform(y_pred_train.reshape(-1, 1))
y_test = y_scaler.inverse_transform(y_test)
y_train = y_scaler.inverse_transform(y_train)

# Evaluation
test_error = np.sqrt(mean_squared_error(y_test, y_pred_test))
train_error = np.sqrt(mean_squared_error(y_train, y_pred_train))

print("Test RMSE:", test_error)
print("Train RMSE:", train_error)
# Calculate absolute percentage error
abs_percentage_error_test = np.abs((y_test - y_pred_test) / y_test)*100

# Calculate mean absolute percentage error
mape_Power = np.mean(abs_percentage_error_test)

print("Mean Absolute Percentage Error (Test):", mape_Power)


# Parity plot for training and test data
plot_limit_min = (df_testing["m_dot_inj"]).min()-2
plot_limit_max = (df_testing["m_dot_inj"]).max()+2
x = np.linspace(plot_limit_min, plot_limit_max, 10)
y = x
y_plus_10 = y * 1.1
y_minus_10 = y * 0.9

plt.scatter(y_train, y_pred_train, alpha=0.9, edgecolor='darkgreen', facecolor='none', marker='^', s=100, linewidths=3, label='Training data')
plt.scatter(y_test, y_pred_test, alpha=0.9, edgecolor='darkblue', facecolor='none', marker='D', s=80, linewidths=2, label='Testing data')

plt.minorticks_on()
plt.tick_params(direction='in', which='minor',bottom=True, top=True, left=True, right=True, length= 3)
plt.tick_params(direction='in', which='major',bottom=True, top=True, left=True, right=True, length= 7)

#plt.colorbar(label="Relative Percent Error",orientation="vertical")
plt.plot(x, y, color='r', alpha=0.4)
plt.plot(x, y_plus_10, linestyle=':', color='k', alpha=0.9, label= " $\pm$ 10 %")
plt.plot(x, y_minus_10, linestyle=':', color='k', alpha=0.9)
plt.xlim(plot_limit_min, plot_limit_max)
plt.ylim(plot_limit_min, plot_limit_max)

plt.xlabel("Actual Injection Mass Flow Rate [lb/min]", fontdict=font)
plt.ylabel("Predicted Injection Mass Flow Rate [lb/min]", fontdict=font)
plt.title("Injection Mass Flow Rate [lb/min]", fontdict=font)

plt.fill_between(x, y, y_plus_10, color='blue', alpha=0.1)

plt.fill_between(x, y, y_minus_10, color='blue', alpha=0.1)


legend_elements = [Line2D([0], [0], marker='D', color='darkblue', linestyle='',
                                          label='Training Data',
                                         markerfacecolor='none', markersize=10),
                   Line2D([0], [0], marker='^', color='darkgreen', linestyle='',
                          label='Testing Data',
                          markerfacecolor='none', markersize=10),
                   Line2D([0], [0], marker='', color='darkgreen', linestyle='',
                          label='MAPE = {:.2f} %'.format(mape_Power),
                          markerfacecolor='none', markersize=10)]
second_legend = plt.legend(handles=legend_elements,loc='upper left', prop={'family': 'Times New Roman', 'size': 12})


plt.tight_layout() 
plt.savefig('Injection mass flow rate', dpi=600)
plt.show() 


#%%


"Discharge Temperature"


plt.figure(figsize=(6,4))
font = {'family': 'Times New Roman', 'color': 'black', 'weight': 'normal', 'size': 14}

# Read data from both sheets
df_training = pd.read_excel('ANN DATA.XLSX', sheet_name=0)
df_testing = pd.read_excel('ANN DATA.XLSX', sheet_name=1)

# Combine data from both sheets
df_combined = pd.concat([df_training, df_testing])

# Split combined data into features and labels
features = df_combined[['P_suc', 'P_dis', 'P_inj', 'F']].to_numpy()

"""
F = compressor Speed 
P_suc = Suction Pressure
P_inj = Injection pressure
P_dis = Discharge Pressure
T_dis = Discharge Temperature"""

labels = df_combined[['T_dis']].to_numpy()

# Scale features and labels
x_scaler = preprocessing.MinMaxScaler()
y_scaler = preprocessing.MinMaxScaler()

X = x_scaler.fit_transform(features)
y = y_scaler.fit_transform(labels)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.8, random_state=1)

# Train the model
neural_network = MLPRegressor(random_state=1, max_iter=3000, hidden_layer_sizes=(20,), activation='relu', solver='lbfgs')
neural_network.fit(X_train, y_train)

# Predictions
y_pred_test = neural_network.predict(X_test)
y_pred_train = neural_network.predict(X_train)

# Inverse scaling for predictions
y_pred_test = y_scaler.inverse_transform(y_pred_test.reshape(-1, 1))
y_pred_train = y_scaler.inverse_transform(y_pred_train.reshape(-1, 1))
y_test = y_scaler.inverse_transform(y_test)
y_train = y_scaler.inverse_transform(y_train)

# Evaluation
test_error = np.sqrt(mean_squared_error(y_test, y_pred_test))
train_error = np.sqrt(mean_squared_error(y_train, y_pred_train))

print("Test RMSE:", test_error)
print("Train RMSE:", train_error)
# Calculate absolute percentage error
abs_percentage_error_test = np.abs((y_test - y_pred_test))

# Calculate mean absolute percentage error
mape_Power = np.mean(abs_percentage_error_test)

print("Mean Absolute Percentage Error (Test):", mape_Power)


# Parity plot for training and test data
plot_limit_min = (df_testing["T_dis"]).min()-2
plot_limit_max = (df_testing["T_dis"]).max()+2
x = np.linspace(plot_limit_min, plot_limit_max, 10)
y = x
y_plus_10 = y * 1.1
y_minus_10 = y * 0.9

plt.scatter(y_train, y_pred_train, alpha=0.9, edgecolor='darkgreen', facecolor='none', marker='^', s=100, linewidths=3, label='Training data')
plt.scatter(y_test, y_pred_test, alpha=0.9, edgecolor='darkblue', facecolor='none', marker='D', s=80, linewidths=2, label='Testing data')

plt.minorticks_on()
plt.tick_params(direction='in', which='minor',bottom=True, top=True, left=True, right=True, length= 3)
plt.tick_params(direction='in', which='major',bottom=True, top=True, left=True, right=True, length= 7)


plt.plot(x, y, color='r', alpha=0.4)
plt.plot(x, y_plus_10, linestyle=':', color='k', alpha=0.9, label= " $\pm$ 10 %")
plt.plot(x, y_minus_10, linestyle=':', color='k', alpha=0.9)
plt.xlim(plot_limit_min, plot_limit_max)
plt.ylim(plot_limit_min, plot_limit_max)

plt.xlabel("Actual Discharge Temperature [K]", fontdict=font)
plt.ylabel("Predicted Discharge Temperature [K]", fontdict=font)
plt.title("Discharge Temperature [K]", fontdict=font)

plt.fill_between(x, y, y_plus_10, color='blue', alpha=0.1)

plt.fill_between(x, y, y_minus_10, color='blue', alpha=0.1)

legend_elements = [Line2D([0], [0], marker='D', color='darkblue', linestyle='',
                                          label='Training Data',
                                         markerfacecolor='none', markersize=10),
                   Line2D([0], [0], marker='^', color='darkgreen', linestyle='',
                          label='Testing Data',
                          markerfacecolor='none', markersize=10),
                   Line2D([0], [0], marker='', color='darkgreen', linestyle='',
                          label='MAE = {:.2f} K'.format(mape_Power),
                          markerfacecolor='none', markersize=10)]
second_legend = plt.legend(handles=legend_elements,loc='upper left', prop={'family': 'Times New Roman', 'size': 12})


plt.tight_layout() 
plt.savefig('Discharge Temperature', dpi=600)
plt.show() 






















