
from pathlib import Path
import csv
import sys

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from scipy.special import j1

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if not (PROJECT_ROOT / 'config' / 'ismConfig.py').is_file():
    raise FileNotFoundError(f'Cannot locate config/ismConfig.py under {PROJECT_ROOT}')
sys.path.insert(0, str(PROJECT_ROOT))
from config.ismConfig import ismConfig

# Existing project data paths (unchanged).
MY_OUTPUT = Path(r'C:\Users\SHARED\EODP_TER_2021\EODP-TS-ISM\myoutput')
REPORT_OUTPUT = Path(r'C:\Users\SHARED\EODP_TER_2021\EODP-TS-ISM\Test_output_ISM_3')

BANDS = ('VNIR-0', 'VNIR-1', 'VNIR-2', 'VNIR-3')
CONTRIBUTORS = ('Diffraction', 'Defocus', 'WFE', 'Detector', 'Smearing', 'Motion')
COLORS = {
    'Diffraction': 'tab:blue', 'Defocus': 'tab:orange',
    'WFE': 'tab:green', 'Detector': 'tab:red',
    'Smearing': 'tab:purple', 'Motion': 'tab:brown',
    'System': 'black',
}


def image_shape(band):
    filepath = MY_OUTPUT / f'ism_toa_optical_{band}.nc'
    if not filepath.is_file():
        raise FileNotFoundError(f'Optical TOA not found: {filepath}')
    with Dataset(filepath, 'r') as ds:
        shape = ds.variables['toa'].shape
    if len(shape) != 2:
        raise ValueError(f'Expected 2D TOA in {filepath}; found {shape}')
    return shape


def calculate_mtf(shape, wavelength, cfg):
    n_alt, n_act = shape
    w = float(cfg.pix_size)
    fn_alt = np.fft.fftshift(np.fft.fftfreq(n_alt))
    fn_act = np.fft.fftshift(np.fft.fftfreq(n_act))
    f_alt, f_act = np.meshgrid(fn_alt, fn_act, indexing='ij')
    fn_2d = np.hypot(f_alt, f_act)
    cutoff = float(cfg.D) / (wavelength * float(cfg.f))
    fr = fn_2d / (w * cutoff)

    hd = np.zeros_like(fr)
    inside = fr < 1.0
    q = fr[inside]
    hd[inside] = (2 / np.pi) * (np.arccos(q) - q * np.sqrt(1 - q ** 2))

    x = np.pi * float(cfg.defocus) * fr * (1 - fr)
    hf = np.ones_like(x)
    valid = np.abs(x) > 1e-12
    hf[valid] = 2 * j1(x[valid]) / x[valid]

    k_wfe = (float(cfg.kLF) * (float(cfg.wLF) / wavelength) ** 2
             + float(cfg.kHF) * (float(cfg.wHF) / wavelength) ** 2)
    hw = np.exp(-fr * (1 - fr) * k_wfe)
    hdet = np.abs(np.sinc(fn_2d))
    hsmear = np.broadcast_to(np.sinc(float(cfg.ksmear) * fn_alt)[:, None], shape)
    hmotion = np.sinc(float(cfg.kmotion) * fn_2d)
    hsys = hd * hf * hw * hdet * hsmear * hmotion

    components = dict(zip(CONTRIBUTORS + ('System',),
                          (hd, hf, hw, hdet, hsmear, hmotion, hsys)))
    return fn_alt, fn_act, components, cutoff


def save_plot(frequencies, curves, band, direction, central_pixel):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    try:
        midpoint = len(frequencies) // 2
        # Plot positive spatial frequencies, including 0.
        x = frequencies[midpoint:]
        for name in CONTRIBUTORS + ('System',):
            ax.plot(x, curves[name][midpoint:],
                    color=COLORS[name],
                    linewidth=2.7 if name == 'System' else 1.4,
                    linestyle='--' if name == 'System' else '-',
                    label=f'{name} MTF')
        ax.axvline(0.5, color='gray', linestyle=':', linewidth=1.4, label='Nyquist')
        ax.set(title=f'System MTF - {direction} - {band} ({central_pixel})',
               xlabel='Normalized spatial frequency [cycles/pixel]',
               ylabel='MTF')
        ax.set_xlim(0, 0.51)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        filename = REPORT_OUTPUT / f'mtf_{direction.lower()}_{band}.png'
        print(f'  Saving: {filename}', flush=True)
        try:
            fig.savefig(filename, dpi=180)
        except OSError as exc:
            # Retain all computed results even if Windows fails to save one image.
            print(f'  WARNING: Could not save {filename.name}: {exc}', flush=True)
            return False
        print(f'  Saved: {filename.name}', flush=True)
        return True
    finally:
        plt.close(fig)


def main():
    print('EODP-ALG-ISM-1030 - System MTF and Nyquist', flush=True)
    print(f'Input folder: {MY_OUTPUT}', flush=True)
    print(f'Output folder: {REPORT_OUTPUT}', flush=True)
    if not MY_OUTPUT.is_dir():
        raise FileNotFoundError(f'MY_OUTPUT not found: {MY_OUTPUT}')
    REPORT_OUTPUT.mkdir(parents=True, exist_ok=True)
    cfg = ismConfig()
    print(f'Pixel size = {cfg.pix_size:.9g} m')
    print(f'D = {cfg.D:.9g} m; f = {cfg.f:.9g} m')
    print('Dominant attenuation is assessed from the smallest individual MTF at Nyquist.')

    rows = []
    failed_plots = []
    for band_index, band in enumerate(BANDS):
        shape = image_shape(band)
        wavelength = float(cfg.wv[band_index])
        fn_alt, fn_act, components, cutoff = calculate_mtf(shape, wavelength, cfg)
        mid_alt, mid_act = shape[0] // 2, shape[1] // 2
        act_curves = {name: mat[mid_alt, :] for name, mat in components.items()}
        alt_curves = {name: mat[:, mid_act] for name, mat in components.items()}

        if not save_plot(fn_act, act_curves, band, 'ACT', f'ALT={mid_alt}'):
            failed_plots.append(f'mtf_act_{band}.png')
        if not save_plot(fn_alt, alt_curves, band, 'ALT', f'ACT={mid_act}'):
            failed_plots.append(f'mtf_alt_{band}.png')

        # For even-sized arrays, Nyquist frequency is -0.5 at array index 0.
        act_values = {name: float(components[name][mid_alt, 0]) for name in CONTRIBUTORS}
        alt_values = {name: float(components[name][0, mid_act]) for name in CONTRIBUTORS}
        defining_act = min(act_values, key=act_values.get)
        defining_alt = min(alt_values, key=alt_values.get)
        nyquist_act = float(components['System'][mid_alt, 0])
        nyquist_alt = float(components['System'][0, mid_act])

        print(f'\n{band} | wavelength={wavelength * 1e9:.1f} nm | shape={shape}')
        print('  Contributor           ACT at Nyquist    ALT at Nyquist')
        for name in CONTRIBUTORS:
            print(f'  {name:<20s} {act_values[name]:>13.6f}    {alt_values[name]:>13.6f}')
        print(f'  {"System":<20s} {nyquist_act:>13.6f}    {nyquist_alt:>13.6f}')
        print(f'  Largest attenuation at Nyquist: ACT={defining_act}; ALT={defining_alt}', flush=True)

        row = {
            'band': band, 'wavelength_nm': wavelength * 1e9,
            'pixel_size_m': float(cfg.pix_size), 'cutoff_frequency_inv_m': cutoff,
            'nyquist_frequency_inv_m': 1 / (2 * float(cfg.pix_size)),
            'mtf_nyquist_ACT': nyquist_act, 'mtf_nyquist_ALT': nyquist_alt,
            'dominant_MTF_ACT': defining_act, 'dominant_MTF_ALT': defining_alt,
        }
        for name in CONTRIBUTORS:
            row[f'{name.lower()}_ACT_nyquist'] = act_values[name]
            row[f'{name.lower()}_ALT_nyquist'] = alt_values[name]
        rows.append(row)

    csv_path = REPORT_OUTPUT / 'ism_mtf_nyquist.csv'
    with csv_path.open('w', newline='', encoding='utf-8') as fp:
        writer = csv.DictWriter(fp, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f'\nCSV saved: {csv_path}', flush=True)
    if failed_plots:
        print('Plots that could not be saved:', ', '.join(failed_plots), flush=True)
    else:
        print('All 8 MTF plots saved successfully.', flush=True)


if __name__ == '__main__':
    main()
