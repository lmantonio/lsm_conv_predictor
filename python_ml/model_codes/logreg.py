import os
import argparse
from datetime import datetime
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold
from sklearn.metrics import make_scorer, fbeta_score, confusion_matrix, balanced_accuracy_score
import matplotlib.pyplot as plt
import matplotlib as mpl
import joblib

from plot_utils import set_roc_curve, set_precision_recall_curve
from result_utils import sweep_threshold, set_best_threshold_metrics, set_summary_results

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
input_path = os.path.join(gen_root, 'dataset', 'ml')
output_path = os.path.join(gen_root, 'results')

# ---------------------------
# Buscar arquivos correspondentes
# ---------------------------
files = os.listdir(input_path)
x_train_file = [f for f in files if args.fragment in f and 'x_train' in f][0]
y_train_file = [f for f in files if args.fragment in f and 'y_train' in f][0]
x_test_file = [f for f in files if args.fragment in f and 'x_test' in f][0]
y_test_file = [f for f in files if args.fragment in f and 'y_test' in f][0]

x_train = pd.read_csv(os.path.join(input_path, x_train_file))
y_train = pd.read_csv(os.path.join(input_path, y_train_file)).squeeze()
x_test = pd.read_csv(os.path.join(input_path, x_test_file))
y_test = pd.read_csv(os.path.join(input_path, y_test_file)).squeeze()

# ---------------------------
# Randomized Search - LogisticRegression com F0.5
# ---------------------------
param_dist = {
    'penalty': ['l2'],
    'C': np.logspace(-4, 4, 20)
}

clf = LogisticRegression(solver='lbfgs', max_iter=5000, random_state=args.random_state)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=args.random_state)

# Scorer F0.5
f0_5_scorer = make_scorer(fbeta_score, beta=0.5)

search = RandomizedSearchCV(
    clf,
    param_distributions=param_dist,
    n_iter=20,
    scoring=f0_5_scorer,
    n_jobs=-1,
    cv=cv,
    verbose=2,
    random_state=args.random_state
)

search.fit(x_train, y_train)

# ---------------------------
# Criar pasta de resultados
# ---------------------------
timestamp = datetime.now().strftime("%Y_%m_%d_%H%M%S")
output_dir_name = f'logreg_{args.input_gen}_{args.fragment}_{timestamp}'
output_dir = os.path.join(output_path, output_dir_name)
general_dir = os.path.join(output_dir, 'general_results')
best_dir = os.path.join(output_dir, 'best_result')
os.makedirs(general_dir, exist_ok=True)
os.makedirs(best_dir, exist_ok=True)

# ---------------------------
# Salvar melhores parâmetros
# ---------------------------
best_model = search.best_estimator_
joblib.dump(best_model, os.path.join(best_dir, 'best_model.pkl'))

best_params_df = pd.DataFrame([search.best_params_])
best_params_df.to_csv(os.path.join(best_dir, 'best_parameters.csv'), index=False)

# ---------------------------
# Avaliar thresholds
# ---------------------------
y_probs = best_model.predict_proba(x_test)[:,1]

results_df = sweep_threshold(y_test, y_probs, beta=0.5)  # <- beta=0.5
results_df.to_csv(os.path.join(general_dir, 'threshold_results.csv'), index=False)

best_thresh_metrics_df = set_best_threshold_metrics(results_df, beta=0.5)
best_thresh_metrics_df.to_csv(os.path.join(best_dir, 'best_threshold_metrics.csv'), index=False)

# ---------------------------
# Salvar curvas
# ---------------------------
set_roc_curve(y_test, y_probs, os.path.join(general_dir, 'roc_curve.png'))
set_precision_recall_curve(y_test, y_probs, os.path.join(general_dir, 'precision_recall_curve.png'))

set_summary_results(output_dir_name, best_thresh_metrics_df, beta=0.5)

print("Processo finalizado")