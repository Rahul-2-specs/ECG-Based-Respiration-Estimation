import wfdb
import numpy as np
import matplotlib.pyplot as plt
import heartpy as hp
from scipy.signal import butter, filtfilt
from scipy.interpolate import interp1d
from scipy.fft import fft, fftfreq

# ============================================================
# ECG-EDR ANALYSIS - MIT-BIH Record 100
# ============================================================

RECORD_NAME = "data/100"
LOWCUT = 0.5
HIGHCUT = 45.0
FILTER_ORDER = 4

ECG_DISPLAY_SECONDS = 10
RESP_MIN_HZ = 0.10
RESP_MAX_HZ = 0.30
RR_RESAMPLE_HZ = 4.0

print("=" * 60)
print("                 ECG-EDR ANALYSIS")
print("=" * 60)

# 1. Load ECG
record = wfdb.rdrecord(RECORD_NAME)
ecg_signal = record.p_signal[:, 0]
fs = float(record.fs)

print(f"Record              : {RECORD_NAME}")
print(f"Signal length       : {len(ecg_signal)} samples")
print(f"Sampling frequency  : {fs:.0f} Hz")
print(f"Duration            : {len(ecg_signal) / fs:.2f} s")


# 2. Band-pass filtering
def butter_bandpass(lowcut, highcut, fs, order=4):
    nyquist = 0.5 * fs
    low = lowcut / nyquist
    high = highcut / nyquist
    return butter(order, [low, high], btype="band")


def bandpass_filter(data, lowcut, highcut, fs, order=4):
    b, a = butter_bandpass(lowcut, highcut, fs, order)
    return filtfilt(b, a, data)


filtered_ecg = bandpass_filter(
    ecg_signal, LOWCUT, HIGHCUT, fs, FILTER_ORDER
)


# 3. R-peak detection
wd, measures = hp.process(filtered_ecg, sample_rate=fs)

r_peaks = np.asarray(wd["peaklist"], dtype=int)
r_peaks = r_peaks[(r_peaks >= 0) & (r_peaks < len(filtered_ecg))]
r_times = r_peaks / fs

if len(r_peaks) < 3:
    raise RuntimeError("Too few R-peaks were detected.")

# 4. RR intervals
rr_intervals = np.diff(r_times)
rr_times = r_times[1:]

mean_rr = np.mean(rr_intervals)
mean_hr = 60.0 / mean_rr

# Keep plausible intervals for interpolation/spectral analysis.
valid = (rr_intervals >= 0.30) & (rr_intervals <= 2.00)
valid_rr_times = rr_times[valid]
valid_rr = rr_intervals[valid]

if len(valid_rr) < 4:
    raise RuntimeError("Too few valid RR intervals for spectral analysis.")


# 5. Uniformly resample RR intervals.
# RR intervals are unevenly sampled, so interpolation is used
# before FFT-based frequency analysis.
uniform_time = np.arange(
    valid_rr_times[0],
    valid_rr_times[-1],
    1.0 / RR_RESAMPLE_HZ
)

if len(uniform_time) < 16:
    raise RuntimeError("RR data are too short for reliable spectral analysis.")

interpolator = interp1d(
    valid_rr_times,
    valid_rr,kind="cubic",
    bounds_error=False,
    fill_value="extrapolate"
)

rr_uniform = interpolator(uniform_time)
rr_detrended = rr_uniform - np.mean(rr_uniform)


# 6. FFT
N = len(rr_detrended)
window = np.hanning(N)
spectrum = np.abs(fft(rr_detrended * window))
freqs = fftfreq(N, 1.0 / RR_RESAMPLE_HZ)

positive = freqs >= 0
freqs = freqs[positive]
spectrum = spectrum[positive]


# 7. Respiration frequency
resp_mask = (freqs >= RESP_MIN_HZ) & (freqs <= RESP_MAX_HZ)

resp_freqs = freqs[resp_mask]
resp_values = spectrum[resp_mask]

if len(resp_freqs) == 0:
    raise RuntimeError("No FFT bins found in the respiration band.")

peak_idx = np.argmax(resp_values)
dominant_frequency = resp_freqs[peak_idx]
respiration_rate_bpm = dominant_frequency * 60.0


# 8. Results
print()
print("-" * 60)
print("                    RESULTS")
print("-" * 60)
print(f"R-peaks detected    : {len(r_peaks)}")
print(f"Mean RR interval    : {mean_rr:.3f} s")
print(f"Mean heart rate     : {mean_hr:.2f} BPM")
print(f"Valid RR intervals  : {len(valid_rr)}")
print(f"Respiration band    : {RESP_MIN_HZ:.2f} - {RESP_MAX_HZ:.2f} Hz")
print(f"Dominant frequency  : {dominant_frequency:.3f} Hz")
print(f"Respiration rate    : {respiration_rate_bpm:.2f} breaths/min")
print("-" * 60)


# ============================================================
# FIGURE 1: Filtered ECG + R-peaks
# ============================================================
display_samples = min(len(filtered_ecg), int(ECG_DISPLAY_SECONDS * fs))
t_ecg = np.arange(display_samples) / fs
display_peaks = r_peaks[r_peaks < display_samples]

plt.figure(figsize=(12, 5))
plt.plot(t_ecg, filtered_ecg[:display_samples],
         linewidth=1.0, label="Filtered ECG")
plt.scatter(display_peaks / fs, filtered_ecg[display_peaks],
            s=28, zorder=3, label="Detected R-peaks")
plt.title("Filtered ECG Signal with R-Peak Detection")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude (mV)")
plt.xlim(0, ECG_DISPLAY_SECONDS)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()


# ============================================================
# FIGURE 2: RR interval variation
# ============================================================
plt.figure(figsize=(12, 5))
plt.plot(rr_times, rr_intervals, marker="o",
         markersize=2.5, linewidth=1.0, label="RR Interval")
plt.title("RR Interval Variation")
plt.xlabel("Time (s)")
plt.ylabel("RR Interval (s)")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()


# ============================================================
# FIGURE 3: Uniformly sampled RR tachogram
# ============================================================
plt.figure(figsize=(12, 5))
plt.plot(uniform_time, rr_uniform,
         linewidth=1.2, label="Interpolated RR Signal")
plt.title("Uniformly Sampled RR Tachogram")
plt.xlabel("Time (s)")
plt.ylabel("RR Interval (s)")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()


# ============================================================
# FIGURE 4: Full low-frequency spectrum
# ============================================================
plt.figure(figsize=(12, 5))
plt.plot(freqs, spectrum, linewidth=1.2, label="RR Spectrum")
plt.axvline(dominant_frequency, linestyle="--", linewidth=1.2,label=f"Dominant = {dominant_frequency:.3f} Hz")
plt.scatter([dominant_frequency], [resp_values[peak_idx]],
            s=45, zorder=3)
plt.title("Respiration Frequency Spectrum from RR Variability")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude")
plt.xlim(0, 0.5)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()


# ============================================================
# FIGURE 5: Respiration band
# ============================================================
plt.figure(figsize=(10, 5))
plt.plot(resp_freqs, resp_values, linewidth=1.5,
         label="Respiration Band")
plt.axvline(dominant_frequency, linestyle="--", linewidth=1.2,
            label=f"{respiration_rate_bpm:.2f} breaths/min")
plt.scatter([dominant_frequency], [resp_values[peak_idx]],
            s=50, zorder=3)
plt.title("Respiration Band (0.10–0.30 Hz)")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude")
plt.xlim(RESP_MIN_HZ, RESP_MAX_HZ)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.show()
