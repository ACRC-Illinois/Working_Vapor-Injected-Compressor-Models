"""
Developed by AMJID KHAN
PhD Student of Dr. Craig Bradshaw
Phsics Inspired Model for Performance Prediction of Vapor Injected Compressor
"""


"Injection mass flow rate model"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error, mean_absolute_percentage_error
from matplotlib.lines import Line2D

#Define compressor non-linear model
def m_dot_inj(x, a0, a1, a2, a3, a4, a5,a6):
    
    rho_inj, F, P_1, P_2, T_inj  = x  
    
    """rho_inj=injection condition density
    F = compressor Speed 
    P_1 = Suction Pressure
    P_2 = Injection pressure
    P_3 = Discharge Pressure
    T_inj = Injection Temperature"""
    
    m_dot_inj = a0+ a1 *(rho_inj) * ((F / F[12]) ** a2) * ((1 - (P_1 / P_2) ** a3)**a4 ) * ((1/T_inj )**a5)
    
    return m_dot_inj

# Load data from Excel
filename = "Import your excel file.xlsx"

df_train = pd.read_excel(filename, sheet_name=0)  # Training set
df_test = pd.read_excel(filename, sheet_name=0)   # Testing set

# Extract independent and dependent variables

x_train = (
    df_train["rho_inj"],
    df_train["F"],
    df_train["P_1"],
    df_train["P_2"],
    df_train["T_inj"]
)
y_train = df_train["m_dot_i"]  # Target variable


x_test = (
    df_test["rho_inj"],
    df_test["F"],
    df_test["P_1"],
    df_test["P_2"],
    df_test["T_inj"]
)
y_test = df_test["m_dot_i"]

# Initial guess for parameters
initial_guess = [0.1, 0.01, 1, 1, 1, 1, 1]  
maxfev = 6000000

# Apply curve fitting to minimize SSE
#For injection mass flow rate, method "trf" works better in most of the cases
params, covariance = curve_fit(m_dot_inj, x_train, y_train, p0=initial_guess, method="lm", maxfev=500000)  

# Extract optimized parameters
a0, a1, a2, a3, a4, a5, a6 = params
print(f"Optimized Parameters: a0={a0}, a1={a1}, a2={a2},  a3={a3}, a4={a4}, a5={a5}, a6={a6}")

# Make predictions on training and test data
y_pred_train = m_dot_inj(x_train, *params)
y_pred_test = m_dot_inj(x_test, *params)

# Compute SSE (Sum of Squared Errors)
SSE_train = np.sum((y_train - y_pred_train) ** 2)
SSE_test = np.sum((y_test - y_pred_test) ** 2)

# Compute R-squared
SS_total_train = np.sum((y_train - np.mean(y_train)) ** 2)
SS_total_test = np.sum((y_test - np.mean(y_test)) ** 2)

R2_train = 1 - (SSE_train / SS_total_train)
R2_test = 1 - (SSE_test / SS_total_test)

rmse = np.sqrt(mean_squared_error(y_test, y_pred_test))
print("RMSE:", rmse)

mse = mean_squared_error(y_test, y_pred_test)
print("MSE", mse)


# Compute MAPE
mape_train = mean_absolute_percentage_error(y_train, y_pred_train) * 100
mape_test = mean_absolute_percentage_error(y_test, y_pred_test) * 100

print(f"SSE (Train): {SSE_train}, SSE (Test): {SSE_test}")
print(f"R² (Train): {R2_train}, R² (Test): {R2_test}")
print(f"MAPE (Train): {mape_train:.2f}%, MAPE (Test): {mape_test:.2f}%")

# Parity plot for training and test data
plot_limit_min = (df_test["m_dot_i"]).min()-1
plot_limit_max = (df_test["m_dot_i"]).max()+1
x = np.linspace(plot_limit_min, plot_limit_max, 10)
y = x
y_plus_10 = y * 1.1
y_minus_10 = y * 0.9

plt.scatter(y_train, y_pred_train, alpha=1, cmap='plasma', facecolor='pink', edgecolors='black', marker='o', s=70, linewidths=1, label='Training data')
plt.scatter(y_test, y_pred_test, alpha=1, cmap='plasma', facecolor='green', 
                        edgecolors='black', marker='o', s=70, linewidths=1, label='Testing data')

plt.minorticks_on()
plt.tick_params(direction='in', which='minor', bottom=True, top=True, left=True, right=True, length=3)
plt.tick_params(direction='in', which='major', bottom=True, top=True, left=True, right=True, length=7)
plt.grid(True, which='major', linestyle='-', linewidth=0.5, color='gray', alpha=0.5)


plt.plot(x, y, color='r', alpha=0.6)
plt.plot(x, y_plus_10, linestyle=':', color='b', alpha=0.6, label="± 10 %")
plt.plot(x, y_minus_10, linestyle=':', color='b', alpha=0.6)
plt.xlim(plot_limit_min, plot_limit_max)
plt.ylim(plot_limit_min, plot_limit_max)

plt.xlabel("Measured Injection Mass Flow Rate")
plt.ylabel("Predicted Injection Mass Flow Rate")
plt.title("Injection Mass Flow Rate")

plt.fill_between(x, y, y_plus_10, color='pink', alpha=0.1)
plt.fill_between(x, y, y_minus_10, color='pink', alpha=0.1)

legend_elements = [
    Line2D([0], [0], marker='o', color='black', linestyle='', label='Training Data', markerfacecolor='pink', markersize=10),
    Line2D([0], [0], marker='o', color='black', linestyle='', label='Testing Data', markerfacecolor='green', markersize=10),
    Line2D([0], [0], marker='_', color='black', linestyle='', label='MAPE = {:.2f} %'.format(mape_test), markerfacecolor='none', markersize=10)
]
plt.legend(handles=legend_elements, loc='upper left', prop={'size': 12})

plt.tight_layout()
plt.savefig('Name your figure as you feel fit', dpi=600)
plt.show()

#%%

"Suction Mass Flow Rate Model"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import mean_absolute_percentage_error,mean_squared_error
from matplotlib.lines import Line2D

# === 1. Define your injection mass flow rate model ===
def m_dot_suc(x, b0, b1, b2):
    rho_suc, F, P3, P1, P2 = x
    
    """rho_suc=Suction condition density
    F = compressor Speed 
    P_1 = Suction Pressure
    P_2 = Injection pressure
    P_3 = Discharge Pressure
    T_inj = Injection Temperature"""
    
    P3_by_P1 = P3 / (P1)
    P2_by_P1 = P2 / (P1)
    
    m_dot_suc = rho_suc * F * (b0 + b1 * P3_by_P1 + b2 * P2_by_P1)
    
    return m_dot_suc

# === 2. Load your training/testing data from Excel ===
filename = "Import your excel file.xlsx"
df_train = pd.read_excel(filename, sheet_name=0)
df_test = pd.read_excel(filename, sheet_name=0)

# === 3. Extract input and output variables ===
x_train = np.vstack((
    df_train["rho_suc"],
    df_train["F"],
    df_train["P_1"],
    df_train["P_2"],
    df_train["P_3"]
))
y_train = df_train["m_dot_suc"]

x_test = np.vstack((
    df_test["rho_suc"],
    df_test["F"],
    df_test["P_1"],
    df_test["P_2"],
    df_test["P_3"]
))
y_test = df_test["m_dot_suc"]

# === 4. Initial parameter guess and curve fitting ===
initial_guess = [0.1, 0.1, 0.1]
params, covariance = curve_fit(m_dot_suc, x_train, y_train, p0=initial_guess, method='lm')

# === 5. Predict on training and testing sets ===
y_pred_train = m_dot_suc(x_train, *params)
y_pred_test = m_dot_suc(x_test, *params)

rmse = np.sqrt(mean_squared_error(y_test, y_pred_test))
print("RMSE:", rmse)
# === 6. Compute error metrics ===
mape_train = mean_absolute_percentage_error(y_train, y_pred_train) * 100
mape_test = mean_absolute_percentage_error(y_test, y_pred_test) * 100

print(f"Optimized Parameters: b0={params[0]:.5f}, b1={params[1]:.5f}, b2={params[2]:.5f}")
print(f"MAPE (Train): {mape_train:.2f}%, MAPE (Test): {mape_test:.2f}%")

# === 7. Parity Plot ===
plot_limit_min = y_test.min() - 2
plot_limit_max = y_test.max() + 5
x = np.linspace(plot_limit_min, plot_limit_max, 100)
y = x
y_plus_10 = y * 1.1
y_minus_10 = y * 0.9

plt.figure(figsize=(6, 4))

plt.scatter(y_train, y_pred_train, alpha=1, cmap='plasma', facecolor='pink', edgecolors='black', marker='o', s=70, linewidths=1, label='Training data')
plt.scatter(y_test, y_pred_test, alpha=1, cmap='plasma', facecolor='green', 
                        edgecolors='black', marker='o', s=70, linewidths=1, label='Testing data')

plt.minorticks_on()
plt.tick_params(direction='in', which='minor', bottom=True, top=True, left=True, right=True, length=3)
plt.tick_params(direction='in', which='major', bottom=True, top=True, left=True, right=True, length=7)
plt.grid(True, which='major', linestyle='-', linewidth=0.5, color='gray', alpha=0.5)


plt.plot(x, y, color='r', alpha=0.6)
plt.plot(x, y_plus_10, linestyle=':', color='b', alpha=0.6, label="± 10 %")
plt.plot(x, y_minus_10, linestyle=':', color='b', alpha=0.6)
plt.xlim(plot_limit_min, plot_limit_max)
plt.ylim(plot_limit_min, plot_limit_max)

plt.xlabel("Measured Suction Mass Flow Rate")
plt.ylabel("Predicted Suction Mass Flow Rate")
plt.title("Suction Mass Flow Rate")

plt.fill_between(x, y, y_plus_10, color='pink', alpha=0.1)
plt.fill_between(x, y, y_minus_10, color='pink', alpha=0.1)

legend_elements = [
    Line2D([0], [0], marker='o', color='black', linestyle='', label='Training Data', markerfacecolor='pink', markersize=10),
    Line2D([0], [0], marker='o', color='black', linestyle='', label='Testing Data', markerfacecolor='green', markersize=10),
    Line2D([0], [0], marker='_', color='black', linestyle='', label='MAPE = {:.2f} %'.format(mape_test), markerfacecolor='none', markersize=10)
]
plt.legend(handles=legend_elements, loc='upper left', prop={'size': 12})

plt.tight_layout()
plt.savefig('Name your figure as you feel fit', dpi=600)
plt.show()





#%%

"Compressor Power Consumption"


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error, mean_absolute_percentage_error
from matplotlib.lines import Line2D

#Define compressor non-linear model
def W_dot_comp(x, c0, c1, c2, c3, c4):
    F, P_1, P_2, P_3 = x  # Ensure correct input unpacking
    
    """
    F = compressor Speed 
    P_1 = Suction Pressure
    P_2 = Injection pressure
    P_3 = Discharge Pressure"""
    
    P3_ref = 1500
    W_dot_comp = c0+ F*(c1*(P_3/P_1)**c4+c2*(P_2/P_1)**c4+c3*(P_3/P3_ref)**c4)
    return  W_dot_comp

# Load data from Excel
filename = "Import your Excel file.xlsx"

df_train = pd.read_excel(filename, sheet_name=0)  # Training set
df_test = pd.read_excel(filename, sheet_name=0)   # Testing set

# Extract independent and dependent variables

x_train = (df_train["F"], df_train["P_1"], df_train["P_2"], df_train["P_3"])  # Features
y_train = df_train["W_dot_comp"]  # Target variable


x_test = (df_test["F"], df_test["P_1"], df_test["P_2"], df_test["P_3"])
y_test = df_test["W_dot_comp"]

# Initial guess for parameters
initial_guess = [0.1, 0.01, 1, 1, 1, 1]  
maxfev = 6000000

# Apply curve fitting to minimize SSE
params, covariance = curve_fit(W_dot_comp, x_train, y_train, p0=initial_guess, method="lm", maxfev=6000000)

# Extract optimized parameters
c0, c1, c2, c3, c4 = params
print(f"Optimized Parameters: c0={c0}, c1={c1}, c2={c2},  c3={c3}, c4={c4}")

# Make predictions on training and test data
y_pred_train = W_dot_comp(x_train, *params)
y_pred_test = W_dot_comp(x_test, *params)

# Compute SSE (Sum of Squared Errors)
SSE_train = np.sum((y_train - y_pred_train) ** 2)
SSE_test = np.sum((y_test - y_pred_test) ** 2)

# Compute R-squared
SS_total_train = np.sum((y_train - np.mean(y_train)) ** 2)
SS_total_test = np.sum((y_test - np.mean(y_test)) ** 2)

R2_train = 1 - (SSE_train / SS_total_train)
R2_test = 1 - (SSE_test / SS_total_test)

# Compute MAPE
mape_train = mean_absolute_percentage_error(y_train, y_pred_train) * 100
mape_test = mean_absolute_percentage_error(y_test, y_pred_test) * 100

rmse = np.sqrt(mean_squared_error(y_test, y_pred_test))
print("RMSE:", rmse)

print(f"SSE (Train): {SSE_train}, SSE (Test): {SSE_test}")
print(f"R² (Train): {R2_train}, R² (Test): {R2_test}")
print(f"MAPE (Train): {mape_train:.2f}%, MAPE (Test): {mape_test:.2f}%")

# Parity plot for training and test data
plot_limit_min = (df_test["W_dot_comp"]).min()-200
plot_limit_max = (df_test["W_dot_comp"]).max()+200
x = np.linspace(plot_limit_min, plot_limit_max, 10)
y = x
y_plus_10 = y * 1.1
y_minus_10 = y * 0.9

plt.scatter(y_train, y_pred_train, alpha=1, cmap='plasma', facecolor='pink', edgecolors='black', marker='o', s=70, linewidths=1, label='Training data')
plt.scatter(y_test, y_pred_test, alpha=1, cmap='plasma', facecolor='green', 
                        edgecolors='black', marker='o', s=70, linewidths=1, label='Testing data')

plt.minorticks_on()
plt.tick_params(direction='in', which='minor', bottom=True, top=True, left=True, right=True, length=3)
plt.tick_params(direction='in', which='major', bottom=True, top=True, left=True, right=True, length=7)
plt.grid(True, which='major', linestyle='-', linewidth=0.5, color='gray', alpha=0.5)


plt.plot(x, y, color='r', alpha=0.6)
plt.plot(x, y_plus_10, linestyle=':', color='b', alpha=0.6, label="± 10 %")
plt.plot(x, y_minus_10, linestyle=':', color='b', alpha=0.6)
plt.xlim(plot_limit_min, plot_limit_max)
plt.ylim(plot_limit_min, plot_limit_max)

plt.xlabel("Measured Power [W]")
plt.ylabel("Predicted Power [W]")
plt.title("Power")

plt.fill_between(x, y, y_plus_10, color='pink', alpha=0.1)
plt.fill_between(x, y, y_minus_10, color='pink', alpha=0.1)

legend_elements = [
    Line2D([0], [0], marker='o', color='black', linestyle='', label='Training Data', markerfacecolor='pink', markersize=10),
    Line2D([0], [0], marker='o', color='black', linestyle='', label='Testing Data', markerfacecolor='green', markersize=10),
    Line2D([0], [0], marker='_', color='black', linestyle='', label='MAPE = {:.2f} %'.format(mape_test), markerfacecolor='none', markersize=10)
]
plt.legend(handles=legend_elements, loc='upper left', prop={'size': 12})

plt.tight_layout()
plt.savefig('Name your figure as you feel fit', dpi=600)
plt.show()


#%%

"Discharge Temperature"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error, mean_absolute_error
from matplotlib.lines import Line2D

#Define compressor non-linear model
def T_dis(x, d0, d1, d2, d3, d4):
    T_s, P_3, P_1, P_2, F = x  # Ensure correct input unpacking
    """T_s=Suction Temperature
    F = compressor Speed 
    P_1 = Suction Pressure
    P_2 = Injection pressure
    P_3 = Discharge Pressure"""    
    
    T_dis = d0+T_s*(d1*(P_3/P_1)**d2+d3*(P_2/P_1)**d2)+d4*F
    
    return T_dis

# Load data from Excel
filename = "Import your excel file.xlsx"

df_train = pd.read_excel(filename, sheet_name=0)  # Training set
df_test = pd.read_excel(filename, sheet_name=0)   # Testing set

# Extract independent and dependent variables

x_train = (
    df_train["T_s"],
    df_train["P_3"],
    df_train["P_1"],
    df_train["P_2"],
    df_train["F"]
)
y_train = df_train["T_dis"]  # Target variable


x_test = (
    df_test["T_s"],
    df_test["P_3"],
    df_test["P_1"],
    df_test["P_2"],
    df_test["F"]
)
y_test = df_test["T_dis"]

# Initial guess for parameters
initial_guess = [0.1, 0.01, 1, 1, 1, 1]  
maxfev = 6000000

# Apply curve fitting to minimize SSE
params, covariance = curve_fit(T_dis, x_train, y_train, p0=initial_guess, method="lm", maxfev=500000)

# Extract optimized parameters
d0, d1, d2, d3, d4 = params
print(f"Optimized Parameters: d0={d0}, d1={d1}, d2={d2}, d3={d3}, d4={d4}")

# Make predictions on training and test data
y_pred_train = T_dis(x_train, *params)
y_pred_test = T_dis(x_test, *params)

# Compute SSE (Sum of Squared Errors)
SSE_train = np.sum((y_train - y_pred_train) ** 2)
SSE_test = np.sum((y_test - y_pred_test) ** 2)

# Compute R-squared
SS_total_train = np.sum((y_train - np.mean(y_train)) ** 2)
SS_total_test = np.sum((y_test - np.mean(y_test)) ** 2)

R2_train = 1 - (SSE_train / SS_total_train)
R2_test = 1 - (SSE_test / SS_total_test)

rmse = np.sqrt(mean_squared_error(y_test, y_pred_test))
print("RMSE:", rmse)


# Compute MAPE
mape_train = mean_absolute_error(y_train, y_pred_train)
mape_test = mean_absolute_error(y_test, y_pred_test)

print(f"SSE (Train): {SSE_train}, SSE (Test): {SSE_test}")
print(f"R² (Train): {R2_train}, R² (Test): {R2_test}")
print(f"MAPE (Train): {mape_train:.2f}%, MAPE (Test): {mape_test:.2f}%")

# Parity plot for training and test data
plot_limit_min = (df_test["T_dis"]).min()-20
plot_limit_max = (df_test["T_dis"]).max()+20
x = np.linspace(plot_limit_min, plot_limit_max, 10)
y = x
y_plus_10 = y * 1.1
y_minus_10 = y * 0.9

plt.scatter(y_train, y_pred_train, alpha=1, cmap='plasma', facecolor='pink', edgecolors='black', marker='o', s=70, linewidths=1, label='Training data')
plt.scatter(y_test, y_pred_test, alpha=1, cmap='plasma', facecolor='green', 
                        edgecolors='black', marker='o', s=70, linewidths=1, label='Testing data')

plt.minorticks_on()
plt.tick_params(direction='in', which='minor', bottom=True, top=True, left=True, right=True, length=3)
plt.tick_params(direction='in', which='major', bottom=True, top=True, left=True, right=True, length=7)
plt.grid(True, which='major', linestyle='-', linewidth=0.5, color='gray', alpha=0.5)


plt.plot(x, y, color='r', alpha=0.6)
plt.plot(x, y_plus_10, linestyle=':', color='b', alpha=0.6, label="± 10 %")
plt.plot(x, y_minus_10, linestyle=':', color='b', alpha=0.6)
plt.xlim(plot_limit_min, plot_limit_max)
plt.ylim(plot_limit_min, plot_limit_max)

plt.xlabel("Measured Discharge Temperature [K]")
plt.ylabel("Predicted Discharge Temperature [K]")
plt.title("Discharge Temperature")

plt.fill_between(x, y, y_plus_10, color='pink', alpha=0.1)
plt.fill_between(x, y, y_minus_10, color='pink', alpha=0.1)

legend_elements = [
    Line2D([0], [0], marker='o', color='black', linestyle='', label='Training Data', markerfacecolor='pink', markersize=10),
    Line2D([0], [0], marker='o', color='black', linestyle='', label='Testing Data', markerfacecolor='green', markersize=10),
    Line2D([0], [0], marker='_', color='black', linestyle='', label='MAE = {:.2f} K'.format(mape_test), markerfacecolor='none', markersize=10)
]
plt.legend(handles=legend_elements, loc='upper left', prop={'size': 12})

plt.tight_layout()
plt.savefig('Name your figure as you feel fit', dpi=600)
plt.show()
