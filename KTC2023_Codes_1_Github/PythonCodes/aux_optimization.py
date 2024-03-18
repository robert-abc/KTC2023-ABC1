import scipy as sp
import numpy as np
import optuna
from PythonCodes import KTCScoring
from PythonCodes import solver

class Objective:
    def __init__(self, inputs_path, groundTruth_path, cnn, sparse_mesh, ref_data):
        self.inputs_path = inputs_path
        self.gts = []
        self.cnn = cnn
        self.sparse_mesh = sparse_mesh
        self.ref_data = ref_data

        for gt_file in groundTruth_path:
            gt = sp.io.loadmat(gt_file)
            gt = gt['truth']
            self.gts.append(gt)

    def __call__(self, trial):
        smooth_lambda = trial.suggest_float("smooth_lambda", 1e-1, 1e2)
        softThersh_perc = trial.suggest_float("softThersh_perc", 0, 0.85)
        cnn_scale = trial.suggest_float("cnn_scale", 0.1, 5)
        morphOpen_radius = trial.suggest_int("morphOpen_radius", 5, 25)

        dic_parameters = {
            'smooth_lambda':smooth_lambda,
            'softThersh_perc':softThersh_perc,
            'cnn_scale':cnn_scale,
            'morphOpen_radius':morphOpen_radius
        }

        score_list = []

        for ilvl in range(1,8):
            for fn in range (len(self.gts)):
                reconstruction = solver.reconstruct(self.inputs_path[fn], ilvl, self.cnn,
                                                     self.sparse_mesh, self.ref_data, dic_parameters)
                
                file_score = KTCScoring.scoringFunction_faster(self.gts[fn], reconstruction)
                score_list.append(file_score)

                trial.report(file_score, fn*(ilvl-1) + fn)

                if trial.should_prune():
                    raise optuna.TrialPruned()

        final_score = np.mean(score_list)

        return final_score