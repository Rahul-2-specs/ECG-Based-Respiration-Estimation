import wfdb
import numpy as np
import matplotlib.pyplot as plt
import heartpy as hp
from scipy.signal import butter, filtfilt
from scipy.fft import fft, fftfreq

# Load record '100' (make sure 100.dat and 100.hea are in the same folder as your Python script)
record = wfdb.rdrecord('100')

# Access ECG signal (use channel 0)
ecg_signal = record.p_signal[:, 0]   # Lead I (or whichever lead is available)

# Sampling frequency
fs = record.fs

print(f"Signal length: {len(ecg_signal)} samples")
print(f"Sampling frequency: {fs} Hz")

# Bandpass filter design
def butter_bandpass(lowcut, highcut, fs, order=4):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def bandpass_filter(data, lowcut=0.5, highcut=45.0, fs=320, order=4):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y

# Apply filter
filtered_ecg = bandpass_filter(ecg_signal, lowcut=0.5, highcut=45.0, fs=fs)

# Plot the filtered ECG waveform
plt.figure(figsize=(12, 4))
t = np.arange(len(filtered_ecg)) / fs  # time axis in seconds
plt.plot(t, filtered_ecg, label='Filtered ECG')
plt.title('Filtered ECG Signal')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude (mV)')
plt.xlim(0, 10)  # Show only the first 10 seconds for clarity
plt.legend()
plt.grid(True)
plt.show()

# Detect R-peaks using heartpy
wd, m = hp.process(filtered_ecg, sample_rate=fs)

# Get R-peak locations
r_peaks = np.array(wd['peaklist'])
r_peaks = r_peaks[r_peaks > 0]  # Remove invalid peaks

# R-peak times (seconds)
r_times = r_peaks / fs

# RR intervals
rr_intervals = np.diff(r_times)

# FFT of RR intervals
N = len(rr_intervals)
T = np.mean(rr_intervals)   # Average time difference
yf = fft(rr_intervals - np.mean(rr_intervals))  # remove DC component
xf = fftfreq(N, T)

# Focus only on positive frequencies (cut off the DC component)
xf = xf[:N//2]  # Positive frequencies
yf = 2.0/N * np.abs(yf[:N//2])  # Double the amplitude to correct for one-sided FFT

# Find the peak frequency in the respiration band (0.1 to 0.3 Hz)
respiration_band_freqs = xf[(xf >= 0.1) & (xf <= 0.3)]  # Focus on 0.1 - 0.3 Hz for respiration
respiration_band_vals = yf[(xf >= 0.1) & (xf <= 0.3)]

# Find the dominant frequency in the respiration band
dominant_frequency = respiration_band_freqs[np.argmax(respiration_band_vals)]

# Convert to breaths per minute (BPM)
respiration_rate_bpm = dominant_frequency * 60

# Print the respiration rate
print(f"Dominant Respiration Frequency: {dominant_frequency:.2f} Hz")
print(f"Respiration Rate: {respiration_rate_bpm:.2f} BPM")

# Plot respiration band (0 to 0.5 Hz)
plt.figure(figsize=(12,4))
plt.plot(xf, yf)
plt.title('Respiration Frequency Spectrum')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.xlim(0, 0.5)  # Human breathing range is ~0.1 to 0.3 Hz
plt.grid(True)
plt.show()
