import re
import os
import argparse
import joblib
import scipy as sp
import tensorflow as tf
import optuna
from PythonCodes.aux_optimization import Objective

# Get input arguments
parser = argparse.ArgumentParser(description=
            'Optimize hyperparameters of the algorithm.')

parser.add_argument('input_path', type=str,
                    help='Path with meausrements to be reconstructed.')
parser.add_argument('groundTruth_path', type=str,
                    help='Path with ground truth images.')
parser.add_argument('study_file', type=str,
                    help='Name of the file to save the optimization study file.')
parser.add_argument('iterations', type=int,
                    help='Number of iterations to run the optimization.')

args = parser.parse_args()

# Get image names
gt_names = sorted(os.listdir(args.groundTruth_path))
r = re.compile(".*\.mat")
gt_names = list(filter(r.match, gt_names))

img_names = sorted(os.listdir(args.input_path))
r2 = re.compile("^(?!ref)")
img_names = list(filter(r.match, img_names))
img_names = list(filter(r2.match, img_names))
print(f"{len(img_names)} images were found.")

for ifile in range(len(gt_names)):
   gt_names[ifile] = os.path.join(args.groundTruth_path, gt_names[ifile])
   img_names[ifile] = os.path.join(args.input_path, img_names[ifile])

n_iterations = args.iterations
study_file = args.study_file

sparse_mesh = sp.io.loadmat('PythonCodes/Mesh_sparse.mat')
cnn = tf.keras.models.load_model('ultimate_cnn1.h5')
ref_data = sp.io.loadmat('TrainingData/ref.mat')

if os.path.isfile(study_file):
    study = joblib.load(study_file)
else:
    study = optuna.create_study(direction="maximize",
                            pruner=optuna.pruners.MedianPruner(
                                n_startup_trials=5, n_warmup_steps=3, interval_steps=1
                                )
                            )
    default_params = {"smooth_lambda": 20,
                  "softThersh_perc": 0.14,
                  "cnn_scale": 2.5,
                  "morphOpen_radius": 12
                  }
    
    study.enqueue_trial(default_params)

study.optimize(Objective(img_names, gt_names, cnn, sparse_mesh, ref_data), n_trials=n_iterations)
joblib.dump(study, study_file)