import numpy as np
import scipy as sp
import cv2
import skimage.morphology as skm
from MiscCodes import KTCMeshing
from MiscCodes import KTCRegularization
from MiscCodes import KTCFwd
from MiscCodes import KTCAux
from MiscCodes import KTCScoring

def reconstruct(input_path, categoryNbr, cnn, sparse_mesh, ref_data, dic_parameters=None):
  if dic_parameters is None:
    smooth_lambda = 63.31
    softThersh_perc = 0.271
    cnn_scale = 2.82
    morphOpen_radius = 10
  else:
    smooth_lambda = dic_parameters['smooth_lambda']
    softThersh_perc = dic_parameters['softThersh_perc']
    cnn_scale = dic_parameters['cnn_scale']
    morphOpen_radius = dic_parameters['morphOpen_radius']

  input_mat = sp.io.loadmat(input_path)

  deltareco_pixgrid_smooth = smooth_prior_reconstruction(input_mat, categoryNbr, sparse_mesh, ref_data, smooth_lambda)
  deltareco_pixgrid_soft = soft_thresholding(deltareco_pixgrid_smooth, softThersh_perc)
  deltareco_pixgrid_cnn = cnn_inference(deltareco_pixgrid_soft, cnn, cnn_scale)
  deltareco_pixgrid_preseg = pre_segmentation(deltareco_pixgrid_cnn)
  deltareco_pixgrid_mopen = morphological_filtering(deltareco_pixgrid_preseg, morphOpen_radius)
  deltareco_pixgrid_otsu = otsu_segmentation(deltareco_pixgrid_mopen)
  
  return deltareco_pixgrid_otsu

def smooth_prior_reconstruction(mat_dict2, categoryNbr, sparse_mesh, ref_data, lambd=20):
  Nel = 32  # number of electrodes
  z = (1e-6) * np.ones((Nel, 1))  # contact impedances
  mat_dict = ref_data #load the reference data
  Injref = mat_dict["Injref"] #current injections
  Uelref = mat_dict["Uelref"] #measured voltages from water chamber
  Mpat = mat_dict["Mpat"] #voltage measurement pattern
  vincl = np.ones(((Nel - 1),76), dtype=bool) #which measurements to include in the inversion
  rmind = np.arange(0,2 * (categoryNbr - 1),1) #electrodes whose data is removed

  #remove measurements according to the difficulty level
  for ii in range(0,75):
    for jj in rmind:
      if Injref[jj,ii]:
        vincl[:,ii] = 0
      vincl[jj,:] = 0

  # load premade finite element mesh (made using Gmsh, exported to Matlab and saved into a .mat file)
  mat_dict_mesh = sparse_mesh
  g = mat_dict_mesh['g'] #node coordinates
  H = mat_dict_mesh['H'] #indices of nodes making up the triangular elements
  elfaces = mat_dict_mesh['elfaces'][0].tolist() #indices of nodes making up the boundary electrodes

  #Element structure
  ElementT = mat_dict_mesh['Element']['Topology'].tolist()
  for k in range(len(ElementT)):
    ElementT[k] = ElementT[k][0].flatten()
  ElementE = mat_dict_mesh['ElementE'].tolist() #marks elements which are next to boundary electrodes
  for k in range(len(ElementE)):
    if len(ElementE[k][0]) > 0:
      ElementE[k] = [ElementE[k][0][0][0], ElementE[k][0][0][1:len(ElementE[k][0][0])]]
    else:
      ElementE[k] = []

  #Node structure
  NodeC = mat_dict_mesh['Node']['Coordinate']
  NodeE = mat_dict_mesh['Node']['ElementConnection'] #marks which elements a node belongs to
  nodes = [KTCMeshing.NODE(coord[0].flatten(), []) for coord in NodeC]
  for k in range(NodeC.shape[0]):
    nodes[k].ElementConnection = NodeE[k][0].flatten()
  elements = [KTCMeshing.ELEMENT(ind, []) for ind in ElementT]
  for k in range(len(ElementT)):
    elements[k].Electrode = ElementE[k]

  #2nd order mesh data
  H2 = mat_dict_mesh['H2']
  g2 = mat_dict_mesh['g2']
  elfaces2 = mat_dict_mesh['elfaces2'][0].tolist()
  ElementT2 = mat_dict_mesh['Element2']['Topology']
  ElementT2 = ElementT2.tolist()
  for k in range(len(ElementT2)):
    ElementT2[k] = ElementT2[k][0].flatten()
  ElementE2 = mat_dict_mesh['Element2E']
  ElementE2 = ElementE2.tolist()
  for k in range(len(ElementE2)):
    if len(ElementE2[k][0]) > 0:
      ElementE2[k] = [ElementE2[k][0][0][0], ElementE2[k][0][0][1:len(ElementE2[k][0][0])]]
    else:
      ElementE2[k] = []

  NodeC2 = mat_dict_mesh['Node2']['Coordinate']  # ok
  NodeE2 = mat_dict_mesh['Node2']['ElementConnection']  # ok
  nodes2 = [KTCMeshing.NODE(coord[0].flatten(), []) for coord in NodeC2]
  for k in range(NodeC2.shape[0]):
    nodes2[k].ElementConnection = NodeE2[k][0].flatten()
  elements2 = [KTCMeshing.ELEMENT(ind, []) for ind in ElementT2]
  for k in range(len(ElementT2)):
    elements2[k].Electrode = ElementE2[k]

  Mesh = KTCMeshing.Mesh(H,g,elfaces,nodes,elements)
  Mesh2 = KTCMeshing.Mesh(H2,g2,elfaces2,nodes2,elements2)

  sigma0 = np.ones((len(Mesh.g), 1)) #linearization point
  corrlength = 1 * 0.115 #used in the prior
  var_sigma = 0.05 ** 2 #prior variance
  mean_sigma = sigma0
  smprior = KTCRegularization.SMPrior(Mesh.g, corrlength, var_sigma, mean_sigma)

  # set up the forward solver for inversion
  solver = KTCFwd.EITFEM(Mesh2, Injref, Mpat, vincl)

  vincl = vincl.T.flatten()

  # set up the noise model for inversion
  noise_std1 = 0.05;  # standard deviation for first noise component (relative to each voltage measurement)
  noise_std2 = 0.01;  # standard deviation for second noise component (relative to the largest voltage measurement)
  solver.SetInvGamma(noise_std1, noise_std2, Uelref)

  # Get info from the data
  Inj = mat_dict2["Inj"]
  Uel = mat_dict2["Uel"]
  Mpat = mat_dict2["Mpat"]
  deltaU = Uel - Uelref

  # Reconstruct image
  Usim = solver.SolveForward(sigma0, z) #forward solution at the linearization point
  J = solver.Jacobian(sigma0, z)

  mask = np.array(vincl, bool)
  deltareco = np.linalg.solve(J.T @ solver.InvGamma_n[np.ix_(mask,mask)] @ J + lambd * smprior.L.T @ smprior.L,J.T @ solver.InvGamma_n[np.ix_(mask,mask)] @ deltaU[vincl])

  # interpolate the reconstruction into a pixel image
  deltareco_pixgrid = KTCAux.interpolateRecoToPixGrid(deltareco, Mesh)

  return deltareco_pixgrid

def soft_thresholding(deltareco_pixgrid, perc = 0.14):
  T = np.max(np.abs(deltareco_pixgrid)) * perc
  xf = np.abs(deltareco_pixgrid) - T
  xf = (xf > 0) * xf
  deltareco_pixgrid = np.sign(deltareco_pixgrid) * xf

  return deltareco_pixgrid

def cnn_inference(deltareco_pixgrid, net, prop_scale = 2.5):
  C = prop_scale * np.max(np.abs(deltareco_pixgrid))
  S_max = C
  S_min = -C
  delta_S = S_max - S_min

  S_max_x = np.max(deltareco_pixgrid)
  S_min_x = np.min(deltareco_pixgrid)

  valor_max_proj = (S_max_x - S_min) / delta_S
  valor_min_proj = (S_min_x - S_min) / delta_S

  deltareco_pixgrid = valor_min_proj + (deltareco_pixgrid - S_min_x)/(S_max_x - S_min_x) * (valor_max_proj - valor_min_proj)

  deltareco_pixgrid = cv2.resize(deltareco_pixgrid, [64, 64], interpolation=cv2.INTER_CUBIC)
  deltareco_pixgrid = deltareco_pixgrid.reshape([1,64,64,1])

  output = net.predict(deltareco_pixgrid)
  deltareco_pixgrid = np.squeeze(output).astype('float64')

  deltareco_pixgrid = cv2.resize(deltareco_pixgrid, [256, 256], interpolation=cv2.INTER_CUBIC)

  return deltareco_pixgrid

def pre_segmentation(deltareco_pixgrid):
  c1 = (np.max(deltareco_pixgrid) - 0.5) / 4
  c2 = (np.min(deltareco_pixgrid) - 0.5) / 4

  deltareco_pixgrid[deltareco_pixgrid > (0.5 + c1)] = 1 # Conductive
  deltareco_pixgrid[deltareco_pixgrid < (0.48 + c2)] = 0 # Resistive

  return deltareco_pixgrid

def morphological_filtering(deltareco_pixgrid, disk_radius = 12):
  se = skm.disk(disk_radius)
  deltareco_pixgrid = skm.opening(deltareco_pixgrid, footprint = se)

  return deltareco_pixgrid

def otsu_segmentation(deltareco_pixgrid):
  # threshold the image histogram using Otsu's method
  level, x = KTCScoring.Otsu2(deltareco_pixgrid.flatten(), 256, 7)

  deltareco_pixgrid_segmented = np.zeros_like(deltareco_pixgrid)

  ind0 = deltareco_pixgrid < x[level[0]]
  ind1 = np.logical_and(deltareco_pixgrid >= x[level[0]],deltareco_pixgrid <= x[level[1]])
  ind2 = deltareco_pixgrid > x[level[1]]
  inds = [np.count_nonzero(ind0),np.count_nonzero(ind1),np.count_nonzero(ind2)]
  bgclass = inds.index(max(inds)) #background class

  match bgclass:
    case 0:
      deltareco_pixgrid_segmented[ind1] = 2
      deltareco_pixgrid_segmented[ind2] = 2
    case 1:
      deltareco_pixgrid_segmented[ind0] = 1
      deltareco_pixgrid_segmented[ind2] = 2
    case 2:
      deltareco_pixgrid_segmented[ind0] = 1
      deltareco_pixgrid_segmented[ind1] = 1

  return deltareco_pixgrid_segmented