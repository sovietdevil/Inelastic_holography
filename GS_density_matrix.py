import numpy as np
import scipy
from abtem.core.energy import energy2wavelength
from tqdm import tqdm

def intensity_profile(dm1d):
    return np.real(np.diag(dm1d))

def propagation_dm1d(dm1d, distance, energy, sampling):
    n_row = dm1d.shape[0]
    wavelength = energy2wavelength(energy)
    k = np.fft.fftfreq(n_row, sampling)
    kp, kpp = np.meshgrid(k, k)
    kp2 = kp**2
    kpp2 = kpp**2
    kdm1d = scipy.fft.fft2(dm1d)
    kdm1d_prop = kdm1d*np.exp(1j*np.pi*wavelength*distance*(kpp2-kp2))
    dm1d_prop = np.fft.ifft2(kdm1d_prop)
    return dm1d_prop

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

def profile_alignment(target, reference, rx, sampling, iter_time=5):
    result1, peak_index1 = phase_corr(reference, target, sampling)
    shifted = drift(target, -rx[peak_index1], sampling)
    for i in range(iter_time):
        result2, peak_index2 = phase_corr(reference, shifted, sampling)
        shifted = drift(shifted, -rx[peak_index2], sampling)
    deltax = rx[peak_index1] - rx[peak_index2]
    return shifted, result2, deltax

def construct_DM1d(rx, stack_1d, defocus, sampling, energy, n_iter=100, drift_list=None, drift_corr=False, check_align=1, align_times=10, stop_align=None):
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
            dm_prop = propagation_dm1d(dm_1d, -df, energy, sampling)
            dm1d_sum += dm_prop
        dm_mean = dm1d_sum / n_df

        for label, df in enumerate(defocus):
            dm_1d = propagation_dm1d(dm_mean, df, energy, sampling)
            if drift_corr and i_iter % check_align == 0 and i_iter <= stop_align:
                drift_image, result, deltax = profile_alignment(stack_1dp[:,label], np.real(np.diag(dm_1d)), rx, sampling, iter_time=align_times)
                drift_list[label] = deltax
                stack_1dp[:,label] = drift_image
            for i in range(n_row):
                dm_1d[i,i] = stack_1dp[i,label]
            dm_1d_rec[label,:,:] = dm_1d

    return dm_mean, drift_list, stack_1dp