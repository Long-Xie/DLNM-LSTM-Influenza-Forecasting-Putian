#############Influenza A
# -*- coding: utf-8 -*-
import psutil
import os
process = psutil.Process(os.getpid())
print(process.memory_info().rss / 1024**2, "MB")


import sys
import sysconfig
print("site-packages :", sysconfig.get_paths()["purelib"])
import numpy as np
import pandas as pd
import tensorflow as tf
print(tf.executing_eagerly())
from sklearn.preprocessing import MinMaxScaler


df = pd.read_excel(r'./df_merge3_joined.xlsx')
df=df.loc[:,['Date','typeA_rate','temp','tempdiff','humidity','precip','windspeed','winddir','solarradiation','uvindex', "DOW", "weekly_detection", "facial_coverings","outsiders_proportion"]]
df['Date'] = pd.to_datetime(df['Date'], format='%Y-%m-%d')
df = df.sort_values('Date')
df.index = df['Date'] 
df_train = df[df['Date'].dt.year < 2023]
df_test = df[df['Date'].dt.year > 2022]
del df_train['Date']
del df_test['Date']
scaler = MinMaxScaler(feature_range=(0,1))
df_for_training = scaler.fit_transform(df_train)
df_for_testing = scaler.transform(df_test)

def create_data(df, n_future, n_past):
    n_feature = df.shape[1]
    x_data, y_data = [], []
    for i in range(n_past, len(df) - n_future + 1):
        x_data.append(np.array(df[i - n_past:i, 0:n_feature]))
        y_data.append(np.array(df[i, 0])) 
    return np.array(x_data), np.array(y_data)
n_future = 1 
n_past = 15
X_train, y_train = create_data(df_for_training, n_future=n_future, n_past=n_past)
X_test, y_test = create_data(df_for_testing, n_future=n_future, n_past=n_past)
seed = 123
np.random.seed(seed)
tf.random.set_seed(seed)

from hyperopt import fmin, tpe, hp, Trials, STATUS_OK
from tensorflow import keras
from tensorflow.keras.layers import LSTM   #GRU
from tensorflow.keras.layers import Dense
from tensorflow.keras.models import Sequential
def gb_mse_cv(params):
    params = { 'units': int(params['units']),
               'activation': params['activation'] , 
               'batch_size': int(params['batch_size']), 
               'learning_rate': params['learning_rate'],
               'epochs': int(params['epochs']),
               'choice':params['choice']
               }
    model = Sequential()
    if params['choice']['layers']== 'one':
       model.add(LSTM(units=params['units'], activation=params['activation'],input_shape=(X_train.shape[1], X_train.shape[2])))
    if params['choice']['layers']== 'two':
        model.add(LSTM(units=params['units'], activation=params['activation'], return_sequences = True ,input_shape=(X_train.shape[1], X_train.shape[2])))
        model.add(LSTM(units=params['choice']['units2'], activation=params['activation'])) 
    if params['choice']['layers']== 'three':
        model.add(LSTM(units=params['units'], activation=params['activation'], return_sequences = True ,input_shape=(X_train.shape[1], X_train.shape[2])))
        model.add(LSTM(units=params['choice']['units22'], return_sequences = True , activation=params['activation'])) 
        model.add(LSTM(units=params['choice']['units3'], activation=params['activation'])) 
    model.add(Dense(1))
    opt = keras.optimizers.Adam(learning_rate = params['learning_rate'])
    model.compile(loss = 'mse',optimizer = opt, metrics=['mse'] )
    model.fit(X_train, y_train, batch_size=params['batch_size'], epochs = params['epochs'], validation_split=0.1)
    score = model.evaluate(X_test, y_test, batch_size=params['batch_size'])
    
    return {'loss': score[0], 'status': STATUS_OK}
space = {   'units': hp.choice('units', [8, 16, 32, 64, 128, 256, 512]),
            'activation':hp.choice('activation',['sigmoid','relu','tanh','swish','elu','softplus']),
            'batch_size' : hp.choice('batch_size', [8, 16, 32, 64, 128, 256, 512]),
            'learning_rate': hp.choice("learning_rate",[0.1,0.01,0.001,0.0001,0.00001,0.000001]),
            'epochs': hp.quniform('epochs', 300, 1000, 1),
            'choice':hp.choice('num_layers', [
                    {'layers':'one', },
                    {'layers':'two',
                    'units2': hp.choice('units2', [8, 16, 32, 64, 128, 256, 512]), 
                    },
                    {'layers':'three',
                     'units22': hp.choice('units22', [8, 16, 32, 64, 128, 256, 512]), 
                     'units3': hp.choice('units3', [8, 16, 32, 64, 128, 256, 512]), 
                    }           ]  ) 
                    }
trials = Trials()
seed = seed
best=fmin(fn=gb_mse_cv,
          space=space, 
          algo=tpe.suggest, 
          max_evals=150, 
          trials=trials,
          rstate = np.random.default_rng(seed))
best


best = {'activation': 4,
        'batch_size': 3,
        'epochs': 880.0,
        'learning_rate': 4,
        'num_layers': 0,
        'units': 4,
        'units22': 0,
        'units3': 0}


np.random.seed(seed)
tf.random.set_seed(seed)



model = Sequential()


if best['num_layers'] == 0:
    model.add(LSTM(
        units=[8, 16, 32, 64, 128, 256, 512][best['units']],
        activation=['sigmoid','relu','tanh','swish','elu','softplus'][best['activation']],
        return_sequences=False,  
        input_shape=(X_train.shape[1], X_train.shape[2])
    ))


elif best['num_layers'] == 1:
    model.add(LSTM(
        units=[8, 16, 32, 64, 128, 256, 512][best['units']],
        activation=['sigmoid','relu','tanh','swish','elu','softplus'][best['activation']],
        return_sequences=True,
        input_shape=(X_train.shape[1], X_train.shape[2])
    ))
    model.add(LSTM(
        units=[8, 16, 32, 64, 128, 256, 512][best['units22']],
        activation=['sigmoid','relu','tanh','swish','elu','softplus'][best['activation']],
        return_sequences=False  
    ))


elif best['num_layers'] == 2:
    model.add(LSTM(
        units=[8, 16, 32, 64, 128, 256, 512][best['units']],
        activation=['sigmoid','relu','tanh','swish','elu','softplus'][best['activation']],
        return_sequences=True,
        input_shape=(X_train.shape[1], X_train.shape[2])
    ))
    model.add(LSTM(
        units=[8, 16, 32, 64, 128, 256, 512][best['units22']],
        activation=['sigmoid','relu','tanh','swish','elu','softplus'][best['activation']],
        return_sequences=True
    ))
    model.add(LSTM(
        units=[8, 16, 32, 64, 128, 256, 512][best['units3']],
        activation=['sigmoid','relu','tanh','swish','elu','softplus'][best['activation']],
        return_sequences=False   
    ))


model.add(Dense(1))


###############################




opt = keras.optimizers.Adam(learning_rate = [0.1,0.01,0.001,0.0001,0.00001,0.000001][best['learning_rate']]) 
model.compile(loss = 'mse',optimizer = opt, metrics = ['mse'] )
history = model.fit(X_train, y_train, batch_size = [8, 16, 32, 64, 128, 256, 512][best['batch_size']], validation_split=0.1, epochs=int(best['epochs']), verbose=1)
np.random.seed(seed)
tf.random.set_seed(seed)
model.evaluate(X_test, y_test, batch_size = [8, 16, 32, 64, 128, 256, 512][best['batch_size']])
y_pred = model.predict(X_test)
import time
today = time.strftime("%Y-%m-%d",time.localtime(time.time()))
out_y = df_train.columns[0]
prediction_copies_array = np.repeat(y_pred,X_train.shape[2], axis=-1)
pred = scaler.inverse_transform(np.reshape(prediction_copies_array,(len(y_pred),X_train.shape[2])))[:,0]
original_copies_array = np.repeat(y_test,X_train.shape[2], axis=-1)
original = scaler.inverse_transform(np.reshape(original_copies_array,(len(y_test),X_train.shape[2])))[:,0]

import datetime
import matplotlib.pyplot as plt

pre_date = pd.to_datetime(min(df_test.index), format='%Y-%m-%d') + datetime.timedelta(days=n_past)
predict = pd.DataFrame(pred)
predict["create_time"] = pd.date_range(pre_date,periods=len(predict))
print("Max date of prediction\n", max(predict["create_time"]))
actual = pd.DataFrame(original) 
actual["create_time"] = pd.date_range(pre_date,periods=len(actual))
print("Max date of prediction\n", max(actual["create_time"] ))
plotdt = predict.merge(actual, on="create_time") 
plotdt.columns = ['pred', 'date', 'real']
plotdt['pred_round'] = plotdt['pred'].clip(0, 1)
plotdt['pred_round'] = plotdt['pred_round'].where(plotdt['pred_round'] >= 0, 0)
min(plotdt['pred_round'])

outputpath = ".\\Influenza A\\" + out_y + "_pred_" + today + ".csv"
plotdt.to_csv(outputpath, sep=',', index=False, header=True)

###########################
original = plotdt['real']
pred = plotdt['pred_round']
from sklearn import metrics
test_new=[]
predict_new=[]
for k in range(len(original)):
    if original[k]!=0:
        test_new.append(original[k])
        predict_new.append(pred[k]) 
def MAPE(y_true, y_pred):
    y_true, y_pred = np.array(y_true), np.array(y_pred)
    return np.mean(np.fabs((y_true - y_pred) / y_true))
def smape(y_true, y_pred): return 2.0 * np.mean(np.abs(y_pred - y_true) / (np.abs(y_pred) + np.abs(y_true))) * 1
print("MAE:" , metrics.mean_absolute_error(original, pred)) 
print("MSE:" , metrics.mean_squared_error(original, pred)) 
print("RMSE:" , np.sqrt(metrics.mean_squared_error(original, pred))) 
print("MAPE:" , MAPE(test_new, predict_new))
print("SMAPE:" , smape(original, pred)) 
#loss
loss = history.history['loss']
val_loss = history.history['val_loss']
Epochs = list(range(int(best['epochs'])))
print(max(loss))
print(max(val_loss))
loss_dt = pd.DataFrame({'loss': history.history['loss'],
                        'val_loss': history.history['val_loss'],
                        'epochs': list(range(int(best['epochs'])))    })
loss_dt.to_csv(".\\Influenza A/" + out_y + "_loss_" + today + ".csv", sep=',', index=False, header=True)


#################SHAP1
import shap
import pandas as pd
import numpy as np


explainer = shap.DeepExplainer(model, X_train)
shap_values = explainer.shap_values(X_test)  


n_features = X_test.shape[2]


shap_values_2D = np.array(shap_values).reshape(-1, n_features)


X_test_2D = X_test.reshape(-1, n_features)



features = ['typeA_rate','temp','tempdiff','humidity','precip','windspeed','winddir','solarradiation','uvindex', "DOW", "weekly_detection", "facial_coverings","outsiders_proportion"]
x_test_2d = pd.DataFrame(X_test_2D, columns=features)


shap.summary_plot(shap_values_2D, x_test_2d, plot_type="bar", show=False, color='#F2BB6B')
plt.gcf().savefig(r"./Influenza A/SHAP_fig1.png", dpi=500, bbox_inches='tight')
plt.close()


###########SHAP2
import shap
import matplotlib.pyplot as plt
import os
import pandas as pd


shap.summary_plot(
    shap_values_2D,
    x_test_2d,
    plot_type="bar",
    show=False,
    color='#F2BB6B'
)


fig = plt.gcf()


plt.xticks(fontproperties='Arial', size=8, color='black')
plt.yticks(fontproperties='Arial', size=8, color='black')
plt.xlabel('mean(|SHAP value|)', fontproperties='Arial', fontsize=8, color='black')


save_path = r"./LSTM20260123/Influenza A"
os.makedirs(save_path, exist_ok=True)


fig.savefig(os.path.join(save_path, "SHAP_fig11.png"), bbox_inches='tight', dpi=500)
plt.close(fig)
######

import matplotlib.pyplot as plt
import shap


shap.summary_plot(shap_values_2D, x_test_2d, show=False)


plt.xticks(fontproperties='Arial', size=8, color='black')
plt.yticks(fontproperties='Arial', size=8, color='black')
plt.xlabel('SHAP value', fontproperties='Arial', fontsize=8, color='black')


fig = plt.gcf()
fig.savefig(r"./Influenza A/SHAP_fig2.png",
            bbox_inches='tight', dpi=500)
plt.close(fig)


###########################Influenza B


import psutil
import os
process = psutil.Process(os.getpid())
print(process.memory_info().rss / 1024**2, "MB")


import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.preprocessing import MinMaxScaler


df = pd.read_excel(r'./df_merge3_joined.xlsx')
df=df.loc[:,['Date','typeB_rate','temp','tempdiff','humidity','precip','windspeed','winddir','solarradiation','uvindex', "DOW", "weekly_detection", "facial_coverings","outsiders_proportion"]]


df['Date'] = pd.to_datetime(df['Date'], format='%Y-%m-%d')
df = df.sort_values('Date')
df.index = df['Date'] 
df_train = df[df['Date'].dt.year < 2023]
df_test = df[df['Date'].dt.year > 2022]
del df_train['Date']
del df_test['Date']
scaler = MinMaxScaler(feature_range=(0,1))
df_for_training = scaler.fit_transform(df_train)
df_for_testing = scaler.transform(df_test)

def create_data(df, n_future, n_past):
    n_feature = df.shape[1]
    x_data, y_data = [], []
    for i in range(n_past, len(df) - n_future + 1):
        x_data.append(np.array(df[i - n_past:i, 0:n_feature]))
        y_data.append(np.array(df[i, 0])) 
    return np.array(x_data), np.array(y_data)
n_future = 1 
n_past = 15
X_train, y_train = create_data(df_for_training, n_future=n_future, n_past=n_past)
X_test, y_test = create_data(df_for_testing, n_future=n_future, n_past=n_past)
seed = 123
np.random.seed(seed)
tf.random.set_seed(seed)

from hyperopt import fmin, tpe, hp, Trials, STATUS_OK
from tensorflow import keras
from tensorflow.keras.layers import LSTM   #GRU
from tensorflow.keras.layers import Dense
from tensorflow.keras.models import Sequential
def gb_mse_cv(params):
    params = { 'units': int(params['units']),
               'activation': params['activation'] , 
               'batch_size': int(params['batch_size']), 
               'learning_rate': params['learning_rate'],
               'epochs': int(params['epochs']),
               'choice':params['choice']
               }
    model = Sequential()
    if params['choice']['layers']== 'one':
       model.add(LSTM(units=params['units'], activation=params['activation'],input_shape=(X_train.shape[1], X_train.shape[2])))
    if params['choice']['layers']== 'two':
        model.add(LSTM(units=params['units'], activation=params['activation'], return_sequences = True ,input_shape=(X_train.shape[1], X_train.shape[2])))
        model.add(LSTM(units=params['choice']['units2'], activation=params['activation'])) 
    if params['choice']['layers']== 'three':
        model.add(LSTM(units=params['units'], activation=params['activation'], return_sequences = True ,input_shape=(X_train.shape[1], X_train.shape[2])))
        model.add(LSTM(units=params['choice']['units22'], return_sequences = True , activation=params['activation'])) 
        model.add(LSTM(units=params['choice']['units3'], activation=params['activation'])) 
    model.add(Dense(1))
    opt = keras.optimizers.Adam(learning_rate = params['learning_rate']) 
    model.compile(loss = 'mse',optimizer = opt, metrics=['mse'] )
    model.fit(X_train, y_train, batch_size=params['batch_size'], epochs = params['epochs'], validation_split=0.1)
    score = model.evaluate(X_test, y_test, batch_size=params['batch_size'])
    
    return {'loss': score[0], 'status': STATUS_OK}
space = {   'units': hp.choice('units', [8, 16, 32, 64, 128, 256, 512]),
            'activation':hp.choice('activation',['sigmoid','relu','tanh','swish','elu','softplus']),
            'batch_size' : hp.choice('batch_size', [8, 16, 32, 64, 128, 256, 512]),
            'learning_rate': hp.choice("learning_rate",[0.1,0.01,0.001,0.0001,0.00001,0.000001]),
            'epochs': hp.quniform('epochs', 300, 1000, 1),
            'choice':hp.choice('num_layers', [
                    {'layers':'one', },
                    {'layers':'two',
                    'units2': hp.choice('units2', [8, 16, 32, 64, 128, 256, 512]), 
                    },
                    {'layers':'three',
                     'units22': hp.choice('units22', [8, 16, 32, 64, 128, 256, 512]), 
                    'units3': hp.choice('units3', [8, 16, 32, 64, 128, 256, 512]), 
                    }           ]  ) 
                    }
trials = Trials()
seed = seed
best=fmin(fn=gb_mse_cv, 
          space=space, 
          algo=tpe.suggest, 
          max_evals=150, 
          trials=trials,
          rstate = np.random.default_rng(seed))
best


best = {'activation': 5,
        'batch_size': 0,
        'epochs': 931.0,
        'learning_rate': 3,
        'num_layers': 1,
        'units': 6,
        'units2': 5}



np.random.seed(seed)
tf.random.set_seed(seed)
tf.compat.v1.disable_v2_behavior()

model = Sequential()
if best['num_layers'] == 0:
    model.add(LSTM(units=[8, 16, 32, 64, 128, 256, 512][best['units']], activation=['sigmoid','relu','tanh','swish','elu','softplus'][best['activation']], return_sequences = True ,input_shape=(X_train.shape[1], X_train.shape[2])))
if best['num_layers'] == 1:
    model.add(LSTM(units=[8, 16, 32, 64, 128, 256, 512][best['units']], activation=['sigmoid','relu','tanh','swish','elu','softplus'][best['activation']], return_sequences = True ,input_shape=(X_train.shape[1], X_train.shape[2])))
    model.add(LSTM(units=[8, 16, 32, 64, 128, 256, 512][best['units2']], activation=['sigmoid','relu','tanh','swish','elu','softplus'][best['activation']])) #
if best['num_layers'] == 2:
    model.add(LSTM(units=[8, 16, 32, 64, 128, 256, 512][best['units']], activation=['sigmoid','relu','tanh','swish','elu','softplus'][best['activation']], return_sequences = True ,input_shape=(X_train.shape[1], X_train.shape[2])))
    model.add(LSTM(units=[8, 16, 32, 64, 128, 256, 512][best['units22']], return_sequences = True , activation=['sigmoid','relu','tanh','swish','elu','softplus'][best['activation']])) 
    model.add(LSTM(units=[8, 16, 32, 64, 128, 256, 512][best['units3']], activation=['sigmoid','relu','tanh','swish','elu','softplus'][best['activation']])) 
model.add(Dense(1))
opt = keras.optimizers.Adam(learning_rate = [0.1,0.01,0.001,0.0001,0.00001,0.000001][best['learning_rate']]) 
model.compile(loss = 'mse',optimizer = opt, metrics = ['mse'] )
history = model.fit(X_train, y_train, batch_size = [8, 16, 32, 64, 128, 256, 512][best['batch_size']], validation_split=0.1, epochs=int(best['epochs']), verbose=1)
np.random.seed(seed)
tf.random.set_seed(seed)
model.evaluate(X_test, y_test, batch_size = [8, 16, 32, 64, 128, 256, 512][best['batch_size']])
y_pred = model.predict(X_test)
import time
today = time.strftime("%Y-%m-%d",time.localtime(time.time()))
out_y = df_train.columns[0]
prediction_copies_array = np.repeat(y_pred,X_train.shape[2], axis=-1)
pred = scaler.inverse_transform(np.reshape(prediction_copies_array,(len(y_pred),X_train.shape[2])))[:,0]
original_copies_array = np.repeat(y_test,X_train.shape[2], axis=-1)
original = scaler.inverse_transform(np.reshape(original_copies_array,(len(y_test),X_train.shape[2])))[:,0]

import datetime
import matplotlib.pyplot as plt

pre_date = pd.to_datetime(min(df_test.index), format='%Y-%m-%d') + datetime.timedelta(days=n_past)
predict = pd.DataFrame(pred)
predict["create_time"] = pd.date_range(pre_date,periods=len(predict))
print("predic最大日期\n", max(predict["create_time"]))
actual = pd.DataFrame(original)
actual["create_time"] = pd.date_range(pre_date,periods=len(actual))
print("actual最大日期\n", max(actual["create_time"] ))
plotdt = predict.merge(actual, on="create_time") 
plotdt.columns = ['pred', 'date', 'real']

plotdt['pred_round'] = plotdt['pred'].clip(0, 1)
plotdt['pred_round'] = plotdt['pred_round'].where(plotdt['pred_round'] >= 0, 0)
min(plotdt['pred_round'])

outputpath = ".\\Influenza B\\" + out_y + "_pred_" + today + ".csv"
plotdt.to_csv(outputpath, sep=',', index=False, header=True)
#
original = plotdt['real']
pred = plotdt['pred_round']
from sklearn import metrics
test_new=[]
predict_new=[]
for k in range(len(original)):
    if original[k]!=0:
        test_new.append(original[k])
        predict_new.append(pred[k]) 
def MAPE(y_true, y_pred):
    y_true, y_pred = np.array(y_true), np.array(y_pred)
    return np.mean(np.fabs((y_true - y_pred) / y_true))
def smape(y_true, y_pred): return 2.0 * np.mean(np.abs(y_pred - y_true) / (np.abs(y_pred) + np.abs(y_true))) * 1
print("MAE:" , metrics.mean_absolute_error(original, pred)) 
print("MSE:" , metrics.mean_squared_error(original, pred)) 
print("RMSE:" , np.sqrt(metrics.mean_squared_error(original, pred))) 
print("MAPE:" , MAPE(test_new, predict_new))
print("SMAPE:" , smape(original, pred)) 
#loss
loss = history.history['loss']
val_loss = history.history['val_loss']
Epochs = list(range(int(best['epochs'])))
print(max(loss))
print(max(val_loss))
loss_dt = pd.DataFrame({'loss': history.history['loss'],
                        'val_loss': history.history['val_loss'],
                        'epochs': list(range(int(best['epochs'])))    })
loss_dt.to_csv("./Influenza B/" + out_y + "_loss_" + today + ".csv", sep=',', index=False, header=True)



#####SHAP1

import shap
import pandas as pd
import numpy as np


explainer = shap.DeepExplainer(model, X_train)
shap_values = explainer.shap_values(X_test)  


n_features = X_test.shape[2]


shap_values_2D = np.array(shap_values).reshape(-1, n_features)


X_test_2D = X_test.reshape(-1, n_features)


features = ['typeB_rate','temp','tempdiff','humidity','precip','windspeed','winddir','solarradiation','uvindex', "DOW", "weekly_detection", "facial_coverings","outsiders_proportion"]
x_test_2d = pd.DataFrame(X_test_2D, columns=features)


shap.summary_plot(shap_values_2D, x_test_2d, plot_type="bar", show=False, color='#F2BB6B')
plt.gcf().savefig(r"./Influenza B/SHAP_fig1.png", dpi=500, bbox_inches='tight')
plt.close()


##########################################
######################################SHAP2

import shap
import matplotlib.pyplot as plt
import os
import pandas as pd


shap.summary_plot(
    shap_values_2D,
    x_test_2d,
    plot_type="bar",
    show=False,
    color='#F2BB6B'
)

fig = plt.gcf()


plt.xticks(fontproperties='Arial', size=8, color='black')
plt.yticks(fontproperties='Arial', size=8, color='black')
plt.xlabel('mean(|SHAP value|)', fontproperties='Arial', fontsize=8, color='black')


save_path = r"./Influenza B"
os.makedirs(save_path, exist_ok=True)


fig.savefig(os.path.join(save_path, "SHAP_fig11.png"), bbox_inches='tight', dpi=500)
plt.close(fig)





###############################################


import matplotlib.pyplot as plt
import shap


shap.summary_plot(shap_values_2D, x_test_2d, show=False)


plt.xticks(fontproperties='Arial', size=8, color='black')
plt.yticks(fontproperties='Arial', size=8, color='black')
plt.xlabel('SHAP value', fontproperties='Arial', fontsize=8, color='black')

fig = plt.gcf()
fig.savefig(r"./Influenza B/SHAP_fig2.png",
            bbox_inches='tight', dpi=500)
plt.close(fig)














