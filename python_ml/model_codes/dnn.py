import os
import argparse
from datetime import datetime
import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers, regularizers
import keras_tuner as kt
import time
from sklearn.metrics import fbeta_score

from plot_utils import set_roc_curve, set_precision_recall_curve
from result_utils import sweep_threshold, set_best_threshold_metrics, set_summary_results

# ---------------------------
# Forçar uso exclusivo de CPU
# ---------------------------
tf.config.set_visible_devices([], 'GPU')
tf.config.threading.set_intra_op_parallelism_threads(13)
tf.config.threading.set_inter_op_parallelism_threads(13)

# ---------------------------
# Argumentos do script
# ---------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--input_gen', required=True, help="Diretório base do dataset: 'gen1_repo' ou 'gen2_repo'")
parser.add_argument('--fragment', required=True, help="Fragmento do nome do dataset para identificação")
parser.add_argument('--random_state', type=int, default=42)
args = parser.parse_args()

# ---------------------------
# Montar caminhos de input e output
# ---------------------------
BASE_ROOT = '/Users/lmantonio/Documents/SSD_backup'
gen_root = os.path.join(BASE_ROOT, args.input_gen + '_repo')
input_path = os.path.join(gen_root, 'dataset', 'dl')
output_path = os.path.join(gen_root, 'results')

# ---------------------------
# Buscar arquivos correspondentes
# ---------------------------
files = os.listdir(input_path)
x_train_file = [f for f in files if args.fragment in f and 'x_train' in f][0]
y_train_file = [f for f in files if args.fragment in f and 'y_train' in f][0]
x_val_file = [f for f in files if args.fragment in f and 'x_val' in f][0]
y_val_file = [f for f in files if args.fragment in f and 'y_val' in f][0]
x_test_file = [f for f in files if args.fragment in f and 'x_test' in f][0]
y_test_file = [f for f in files if args.fragment in f and 'y_test' in f][0]

x_train = pd.read_csv(os.path.join(input_path, x_train_file)).astype("float32")
y_train = pd.read_csv(os.path.join(input_path, y_train_file)).squeeze().astype("float32")
x_val = pd.read_csv(os.path.join(input_path, x_val_file)).astype("float32")
y_val = pd.read_csv(os.path.join(input_path, y_val_file)).squeeze().astype("float32")
x_test = pd.read_csv(os.path.join(input_path, x_test_file)).astype("float32")
y_test = pd.read_csv(os.path.join(input_path, y_test_file)).squeeze().astype("float32")

# ---------------------------
# Callback customizado para registrar F0.5 no tuner
# ---------------------------
class F05MetricCallback(tf.keras.callbacks.Callback):
    def __init__(self, validation_data):
        super().__init__()
        self.validation_data = validation_data

    def on_epoch_end(self, epoch, logs=None):
        val_x, val_y = self.validation_data
        y_pred = (self.model.predict(val_x, verbose=0) > 0.5).astype(int)
        f05 = fbeta_score(val_y, y_pred, beta=0.5)
        logs['val_f05_score'] = f05
        print(f" — val_f0.5_score: {f05:.4f}")

f05_callback = F05MetricCallback(validation_data=(x_val.values, y_val.values))

# ---------------------------
# Função para construir modelo
# ---------------------------
def build_model(hp):
    model = keras.Sequential()
    model.add(layers.Input(shape=(x_train.shape[1],)))
    
    for i in range(hp.Int("num_layers", 1, 5)):
        # Regularização L1/L2
        l1 = hp.Float(f"l1_{i}", 0.0, 0.01, step=0.001)
        l2 = hp.Float(f"l2_{i}", 0.0, 0.01, step=0.001)
        
        model.add(layers.Dense(
            units=hp.Choice(f"units_{i}", [16, 32, 64, 128, 256]),
            activation=hp.Choice("activation", ["relu", "elu", "tanh"]),
            kernel_initializer=hp.Choice("kernel_initializer", ["he_normal", "glorot_uniform", "lecun_normal"]),
            kernel_regularizer=regularizers.L1L2(l1=l1, l2=l2)
        ))
        
        if hp.Boolean(f"use_batchnorm_{i}"):
            model.add(layers.BatchNormalization())
        model.add(layers.Dropout(rate=hp.Float(f"dropout_{i}", 0.0, 0.5, step=0.1)))
    
    model.add(layers.Dense(1, activation="sigmoid"))
    
    # Optimizer como hiperparâmetro
    optimizer_choice = hp.Choice("optimizer", ["adam", "rmsprop", "sgd"])
    lr = hp.Choice("learning_rate", [1e-2, 1e-3, 1e-4])
    if optimizer_choice == "adam":
        optimizer = keras.optimizers.Adam(learning_rate=lr)
    elif optimizer_choice == "rmsprop":
        optimizer = keras.optimizers.RMSprop(learning_rate=lr)
    else:
        optimizer = keras.optimizers.SGD(learning_rate=lr)
    
    model.compile(
        optimizer=optimizer,
        loss="binary_crossentropy",
        metrics=[
            keras.metrics.BinaryAccuracy(name="accuracy"),
            keras.metrics.Precision(name="precision"),
            keras.metrics.Recall(name="recall"),
            keras.metrics.AUC(name="auc"),
            keras.metrics.AUC(name="prc", curve="PR")
        ]
    )
    
    return model

# ---------------------------
# Criar pastas do experimento
# ---------------------------
timestamp = datetime.now().strftime("%Y_%m_%d_%H%M%S")
root_dir_name = f'dnn_{args.input_gen}_{args.fragment}_{timestamp}'
root_dir = os.path.join(output_path, root_dir_name)
general_dir = os.path.join(root_dir, "general_results")
best_dir = os.path.join(root_dir, "best_result")
tuner_dir = os.path.join(root_dir, "keras_tuner_dir")

os.makedirs(general_dir, exist_ok=True)
os.makedirs(best_dir, exist_ok=True)
os.makedirs(tuner_dir, exist_ok=True)

# ---------------------------
# Ajuste do batch size como hiperparâmetro
# ---------------------------
hp_global = kt.HyperParameters()
batch_size = hp_global.Choice("batch_size", [32, 64, 128, 256])

train_dataset = tf.data.Dataset.from_tensor_slices((x_train.values, y_train.values))
train_dataset = train_dataset.shuffle(buffer_size=1024).batch(batch_size).prefetch(tf.data.AUTOTUNE)

val_dataset = tf.data.Dataset.from_tensor_slices((x_val.values, y_val.values))
val_dataset = val_dataset.batch(batch_size).prefetch(tf.data.AUTOTUNE)

test_dataset = tf.data.Dataset.from_tensor_slices((x_test.values, y_test.values))
test_dataset = test_dataset.batch(batch_size).prefetch(tf.data.AUTOTUNE)

# ---------------------------
# Configurar tuner Hyperband com F0.5 como objetivo
# ---------------------------
tuner = kt.Hyperband(
    build_model,
    objective=kt.Objective("val_f05_score", direction="max"),
    max_epochs=hp_global.Choice("max_epochs", [20, 40, 60]),
    factor=3,
    directory=tuner_dir,
    project_name=f"dnn_{args.input_gen}_{args.fragment}"
)
tuner.oracle.max_trials = 100

early_stop = keras.callbacks.EarlyStopping(
    monitor="val_loss", patience=5, restore_best_weights=True
)

# ---------------------------
# Rodar busca
# ---------------------------
start_time = time.time()
tuner.search(
    train_dataset,
    validation_data=val_dataset,
    epochs=60,
    callbacks=[early_stop, f05_callback],
    verbose=2
)
elapsed = time.time() - start_time

best_hp = tuner.get_best_hyperparameters(1)[0]
best_model = tuner.get_best_models(1)[0]

# ---------------------------
# Salvar melhores parâmetros
# ---------------------------
params_dict = {hp: best_hp.get(hp) for hp in best_hp.values.keys()}
params_df = pd.DataFrame([params_dict])
params_df.to_csv(os.path.join(best_dir, "best_parameters.csv"), index=False)

# ---------------------------
# Salvar modelo
# ---------------------------
best_model.save(os.path.join(best_dir, "best_model.h5"))

# ---------------------------
# Avaliar thresholds
# ---------------------------
y_probs = best_model.predict(test_dataset).ravel()
results_df = sweep_threshold(y_test, y_probs)
results_df.to_csv(os.path.join(general_dir, "threshold_results.csv"), index=False)

best_thresh_metrics_df = set_best_threshold_metrics(results_df)
best_thresh_metrics_df.to_csv(os.path.join(best_dir, "best_threshold_metrics.csv"), index=False)

# ---------------------------
# Salvar curvas
# ---------------------------
set_roc_curve(y_test, y_probs, os.path.join(general_dir, "roc_curve.png"))
set_precision_recall_curve(y_test, y_probs, os.path.join(general_dir, "precision_recall_curve.png"))

set_summary_results(root_dir_name, best_thresh_metrics_df)

print(f"\n✅ Busca completa finalizada! Tempo total: {elapsed/60:.2f} minutos")