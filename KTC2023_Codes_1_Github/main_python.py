import re
import os
import argparse
import scipy as sp
import tensorflow as tf
from PythonCodes import solver

# Get input arguments
parser = argparse.ArgumentParser(description=
            'Reconstruct EIT images with different amounts of missing data.')

parser.add_argument('input_path', type=str,
                    help='Path with meausrements to be reconstructed.')
parser.add_argument('output_path', type=str,
                    help='Path to save the segmentation of the reconstructed images.')
parser.add_argument('category', type=int,
                    choices=range(1,8), metavar='[1-7]',
                    help='Difficulty category number.')

args = parser.parse_args()

# Get image names
img_names = sorted(os.listdir(args.input_path))
r = re.compile(".*\.mat")
img_names = list(filter(r.match,img_names))
print(f"{len(img_names)} images were found.")

sparse_mesh = sp.io.loadmat('Mesh_sparse.mat')
cnn = tf.keras.models.load_model('ultimate_cnn1.h5')
ref_data = sp.io.loadmat('TrainingData/ref.mat')

for img_name in img_names:
    path_in = os.path.join(args.input_path, img_name)
    path_out = os.path.join(args.output_path, img_name)

    reconstruction = solver.reconstruct(path_in, args.category, cnn, sparse_mesh, ref_data)
    sp.io.savemat(path_out, {'reconstruction':reconstruction})
