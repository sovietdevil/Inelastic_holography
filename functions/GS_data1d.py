import numpy as np
import scipy
from abtem.core.energy import energy2wavelength
from tqdm import tqdm
import torch
from scipy.interpolate import interpn
from matplotlib import pyplot as plt
from scipy.linalg import svd

def intensity_profile(dm1d):
    return np.real(np.diag(dm1d))

class Propagation_layer(torch.nn.Module):
    def __init__(self, init_defocus, sampling, energy):
        super().__init__()
        self.sampling = sampling
        self.wavelength = energy2wavelength(energy)
        #self.defocus = init_defocus
        self.defocus = torch.nn.Parameter(torch.tensor(init_defocus).to(torch.float32), requires_grad=True)
        self.shift = torch.nn.Parameter(torch.zeros(1).to(torch.float32), requires_grad=True)

    def forward(self, density_matrix):
        n_row = density_matrix.shape[0]
        wavelength = self.wavelength
        defocus = self.defocus
        k = np.fft.fftfreq(n_row, self.sampling)
        kp, kpp = np.meshgrid(k, k)
        kp2 = kp**2; kpp2 = kpp**2
        kp = torch.tensor(kp); kpp = torch.tensor(kpp)
        kp2 = torch.tensor(kp2); kpp2 = torch.tensor(kpp2)
        dx = self.shift
        shift_xy = torch.exp(2j*torch.pi*(kp*dx - kpp*dx))
        shift_z = torch.exp(1j*torch.pi*wavelength*defocus*(kp2-kpp2))
        dm_k = torch.fft.fft2(density_matrix)
        return torch.fft.ifft2(dm_k * shift_xy * shift_z)

def propagation_dm1d(dm1d, distance, energy, sampling, apply_tcc=False, alpha=1, Cs=1, delta=1):
    n_row = dm1d.shape[0]
    wavelength = energy2wavelength(energy)
    k = np.fft.fftfreq(n_row, sampling)
    kp, kpp = np.meshgrid(k, k)
    kp2 = kp**2
    kpp2 = kpp**2
    kdm1d = scipy.fft.fft2(dm1d)
    if not apply_tcc:
        kdm1d_prop = kdm1d*np.exp(1j*np.pi*wavelength*distance*(kpp2-kp2))
    else:
        phase_factor = np.exp(1j*np.pi*wavelength*distance*(kpp2-kp2))
        E_delta = np.exp(-1/2*(np.pi*wavelength*delta)**2*(kp**2-kpp**2)**2)
        E_alpha = np.exp(-(np.pi*alpha/wavelength)**2*(Cs*wavelength**3*(kp**3-kpp**3)+distance*wavelength*(kp-kpp))**2)
        kdm1d_prop = kdm1d*phase_factor*E_alpha*E_delta
    dm1d_prop = np.fft.ifft2(kdm1d_prop)
    return dm1d_prop

def focal_spread_results(dm1d, defocus_variation, energy, sampling, focal_step=0.5):
    defocus_range = np.arange(-2*defocus_variation, 2*defocus_variation, focal_step)
    weighting_factor = np.exp(-defocus_range/(2*defocus_variation))
    result = 0
    for index, defocus in enumerate(defocus_range):
        result = result + weighting_factor[index] * intensity_profile(propagation_dm1d(dm1d, defocus, energy, sampling))
    return result

def standardize(profile):
    return (profile - np.mean(profile)) / np.std(profile)

def phase_corr(profile1, profile2, sampling):
    p1_FT = np.fft.fft(standardize(profile1))
    p2_FT = np.fft.fft(standardize(profile2))
    k = np.fft.fftfreq(profile1.shape[0], sampling)
    prod  = p1_FT * p2_FT.conj()
    result = np.fft.ifft(prod/np.abs(prod))
    peak = np.argmax(np.abs(result))
    return result, peak

def drift(y, dx, sampling):
    ft_y = np.fft.fft(y)
    kx = np.fft.fftfreq(len(y), sampling)
    drift_factor = np.exp(2j*np.pi*(kx * dx))
    ft_drift = ft_y * drift_factor
    return np.real(np.fft.ifft(ft_drift))

def profile_drift(target, drift_distance, sampling):
    drifted_profile = []
    for i in range(len(drift_distance)):
        drifted = drift(target[:,i], drift_distance, sampling)
        drifted_profile.append(drifted)
    return np.array(drifted_profile)

def profile_alignment(target, reference, rx, sampling, iter_time=5):
    result1, peak_index1 = phase_corr(reference, target, sampling)
    shifted = drift(target, -rx[peak_index1], sampling)
    for i in range(iter_time):
        result2, peak_index2 = phase_corr(reference, shifted, sampling)
        shifted = drift(shifted, -rx[peak_index2], sampling)
    deltax = rx[peak_index1] - rx[peak_index2]
    return shifted, result2, deltax

def eigenvalues_selection(density_matrix, index=50):
    U, S, V = svd(density_matrix)
    S[index:] = 0
    return U @ np.diag(S) @ V

def construct_DM1d(rx, stack_1d, defocus, sampling, energy, n_iter=100, drift_list=None, drift_corr=False, 
                   check_align=1, align_times=10, stop_align=None, apply_tcc=False, alpha=1, Cs=1, delta=1,
                   filter_component=False, index_sel=50):
    n_row, n_df = stack_1d.shape
    dm_1d_rec = np.zeros((n_df, n_row, n_row)).astype(np.complex128)
    stack_1dp = stack_1d.copy()
    if drift_list is None:
        drift_list = np.zeros((n_df,))
    if stop_align is None:
        stop_align = n_iter

    for label, df in enumerate(defocus):
        image = stack_1dp[:,label]
        dm_1d = np.outer(np.sqrt(image), np.sqrt(image))
        dm_1d_rec[label,:,:] = dm_1d
    
    for i_iter in tqdm(range(n_iter)):
        dm1d_sum = np.zeros((n_row, n_row)).astype(np.complex128)
        for label, df in enumerate(defocus):
            dm_1d = dm_1d_rec[label,:,:]
            dm_prop = propagation_dm1d(dm_1d, -df, energy, sampling, apply_tcc, alpha, Cs, delta)
            dm1d_sum += dm_prop
        dm_mean = dm1d_sum / n_df
        if filter_component:
            dm_mean = eigenvalues_selection(dm_mean, index=index_sel)

        for label, df in enumerate(defocus):
            dm_1d = propagation_dm1d(dm_mean, df, energy, sampling, apply_tcc, alpha, Cs, delta)
            if drift_corr and i_iter % check_align == 0 and i_iter <= stop_align:
                drift_image, result, deltax = profile_alignment(stack_1dp[:,label], np.real(np.diag(dm_1d)), rx, sampling, iter_time=align_times)
                drift_list[label] = deltax
                stack_1dp[:,label] = drift_image
            for i in range(n_row):
                dm_1d[i,i] = stack_1dp[i,label]
            dm_1d_rec[label,:,:] = dm_1d

    return dm_mean, drift_list, stack_1dp

def interp_dm(dm_rec, rx_sim, rx_ex, method='splinef2d'):
    Rx_ex, Rxp_ex = np.meshgrid(rx_ex, rx_ex)
    dm_re = np.real(dm_rec)
    dm_im = np.imag(dm_rec)
    dm_interp_re = interpn((rx_sim, rx_sim), dm_re, (Rx_ex, Rxp_ex), method=method)
    dm_interp_im = interpn((rx_sim, rx_sim), dm_im, (Rx_ex, Rxp_ex), method=method)
    density_matrix = dm_interp_re + 1j*dm_interp_im
    return density_matrix

def standardize(profile):
    return (profile - np.mean(profile))/np.std(profile)

def contrast(orig_profile):
#    profile = (orig_profile - np.mean(orig_profile))/np.std(orig_profile)
    profile = orig_profile
    return (np.max(profile) - np.min(profile))/(np.max(profile) + np.min(profile))

def contrast_profile(defocus, focal_series, normalization=True):
    contrast_prof = []
    for n, df in enumerate(defocus):
        contrast_prof.append(contrast(focal_series[:,n]))
    contrast_prof = np.array(contrast_prof)
    if normalization:
        contrast_prof = contrast_prof/np.max(contrast_prof)
    return contrast_prof

def sim_ex_comparisons(rx, density_matrix, focal_series, defocus, energy, sampling):
    for n_df, df in enumerate(defocus):
        dm_prop = propagation_dm1d(density_matrix, df, energy, sampling)
        line_sim = intensity_profile(dm_prop)
        line_ex = focal_series[:, n_df]
        plt.plot(rx, standardize(line_sim))
        plt.plot(rx, standardize(line_ex))
        plt.title(f"defocus {df/10} nm")
        plt.show()