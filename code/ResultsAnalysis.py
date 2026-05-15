#!/usr/bin/env python3
# Copyright (c) 2026 Paolo Rissone and Federico Ricci-Tersenghi


import matplotlib.pyplot as plt
import matplotlib.pylab as plb
import numpy as np
import math
import sys
import itertools
import os
import re
import pandas as pd
import seaborn as sns
from scipy.stats import gaussian_kde
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator
from io import StringIO
import random
from itertools import groupby
import operator




# ----------------------------------------------------------
# ----------------- GLOBAL PARAMETERS ----------------------

B = 1000 #number of bootstrap resampling
tolerance = 0. #Tolerance of results [0,1] (default = 0.)

#Initialize plots
point_type = ['v', '^', '<', '>', 'd', 'h', '8']
clrBP = [ 'blue', 'purple', 'slategray', 'magenta', 'lime', 'yellow', 'cyan', 'pink', 'gold', 'teal', 'navy', 'indigo', 'chocolate', 'darkcyan', 'olive', 'mediumvioletred', 'deepskyblue', 'turquoise', 'peru', 'orchid', 'steelblue', 'darkkhaki', 'brown', 'black' ]
clrMC = ['red', 'firebrick', 'darkred', 'rebeccapurple', 'darkorchid', 'violet']

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.color'] = plt.rcParams['xtick.color']
plt.rcParams['ytick.color'] = plt.rcParams['ytick.color']

plt.rcParams['figure.dpi'] = 300  # Good balance of quality vs size
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.compression'] = 6  # Maximum PDF compression


# ----------------------------------------------------------
# --------------------- FUNCTIONS --------------------------

#Get overall best energy
def extract_optimal_energies(filenames):
    energies = []
    for file in filenames:
        with open(file, 'r') as file:
            lines = file.readlines()
        for line in reversed(lines):
            if line.strip().startswith('#Optimal configuration has E ='):
                energies.append( float(line.split('=', 1)[1].strip()) )
                break
    return energies
    
    
def import_data(filename):
    data = pd.read_csv(filename, delim_whitespace=True, index_col=None, header=None, comment='#', on_bad_lines='skip')
    
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    # Extract column names from header - check if first line starts with "#Graph"
    header_line = lines[0].lstrip('#').strip()
    if lines[0].startswith('#Graph'):
        header_line = lines[1].lstrip('#').strip()
    
    column_names = [col.strip() for col in re.split(r'\s+', header_line)]
    data.columns = column_names
#    data['path'] = os.path.dirname(filename)
#    data['file'] = os.path.basename(filename)
    
    if "maxIters" in data.columns:
        data.rename(columns={"maxIters": "Iter"}, inplace=True)
        
    return data
    

#Extract algo info
def get_info(filename, data, Egs, N):
    
    ext_to_algo = [
        (('.eo',), None, 'EO'),
        (('.taueo',), None, 'tauEO'),
        (('.sa',), None, 'SA'),
        (('.pa',), (r'R1\b',), 'SA'),   # .pa AND R1 (not R10, R11, etc.)
        (('.pa',), None, 'PA'),
        (('.cea',), None, 'CEA'),
        (('.bps', '.bp'), None, 'BP'),
        (('.bpa',), None, 'BPA'),
        (('.bpr', '.bppa'), None, 'BPR'),
    ]

    algo = next(
        (algo for exts, required, algo in ext_to_algo
         if any(filename.endswith(ext) for ext in exts)
         and (required is None or all(re.search(r, filename) for r in required))),
        ''
    )
    
    # Extract parameters using regex patterns
    param_patterns = [
#        ('T0', r'_T0(\d+\.\d+)'),
#        ('T', r'_T(\d+\.\d+)'),
#        ('R', r'_R(\d+)'),
#        ('dh', r'_dh(\d+\.\d+)'),
        ('betaf', r'_betaf(\d+\.\d+)'),  # Check betaf before beta
        ('beta', r'_beta(\d+\.\d+)'),
    ]
        
    param_name = 'beta' #Default values for BP old code (no input beta required)
    param_value = '0.00'
    for keyword, pattern in param_patterns:
        if keyword in filename:
            # Skip beta if betaf is found
            if keyword == 'beta' and 'betaf' in filename:
                continue
            if algo == 'SA' and keyword == 'R':
                continue
            match = re.search(pattern, filename)
            if match:
                param_name = keyword
                param_value = match.group(1)
                break
                    
    # Extract 'f' parameter
    f = np.nan
    if '.bp' in filename:
        f = re.search(r'_f([^_]+)_', filename).group(1)
        
    # Read footer lines from bottom up until first non-# line
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    footer_lines = []
    for line in reversed(lines):
        if line.strip().startswith('#'):
            footer_lines.append(line)
        elif line.strip() == '':  # Skip blank lines
            continue
        else:
            break
    footer_lines.reverse()  # Put back in original order
    
    # Extract info from footer
    replicas = 1
    time = np.nan
    runs = np.nan
    bestE = np.nan
    avgITERS = np.nan
    accepted = np.nan
    spins = []
    
    for footer in footer_lines:
        if footer.startswith('#Execution time: '):
            time_str = footer.split(':', 1)[1].strip()
            time = float(time_str.split()[0])  # Extract float before 'seconds'
        if footer.startswith('#No. runs: ') or footer.startswith('#Initial runs no.: ') or footer.startswith('#Runs no.: '):
            runs = int(footer.split(':', 1)[1].strip())
        if footer.startswith('#No. replicas per run: '):
            replicas = int(footer.split(':', 1)[1].strip())
        if footer.startswith('#Average iterations per run: '):
            avgITERS = float(footer.split(':', 1)[1].strip())
        if footer.startswith('#Average accepted moves: '):
            accepted = float(footer.split(':', 1)[1].strip())
        if footer.startswith('#Optimal configuration has E ='):
            bestE = float(footer.split('=', 1)[1].strip())
        if footer.startswith('# ') and not any(footer.startswith(prefix) for prefix in ['#Execution', '#No.']):
            spin_values = footer.split(' ', 1)[1].strip().split()
            spins.extend([int(val) for val in spin_values if val.isdigit()])
    
    #Extract stats
    maxIters = data[filename]['Iter'].iloc[-1]
    Results = data[filename][ data[filename]['Iter'] == maxIters ]

    avg_bestE = Results['bestE'].mean()
    std_bestE = Results['bestE'].std()
    avg_bestIter = Results['bestIter'].mean()
    avg_totIters = Results['Iter'].mean()

    #Compute <|S|>, <q>, <f>, and scale factor from EVO files
    avg_f = np.nan
    avg_BPit = np.nan
    if '.bp' in filename:
        filenameEVO = re.sub(r'\.bp', r'EVO.bp', filename)
        avg_f = data[filenameEVO]['F'].mean()
        avg_BPit = data[filenameEVO]['itBP'].mean()
        avg_q = data[filenameEVO]['q'].mean()
        avg_sizeS = data[filenameEVO]['|S|'].mean()
        scale_factor = data[filenameEVO]['itBP'].mean() * avg_sizeS  * replicas
    if filename.endswith(('.sa', '.pa')):
        avg_q = 1.
        avg_sizeS = 1.
        scale_factor = N * replicas
    if filename.endswith(('.eo', '.taueo')):
        avg_q = 1.
        avg_sizeS = 1.
        scale_factor = 1.0
    if '.cea' in filename:
        avg_q = Results['<q>'].mean()
        avg_sizeS = Results['<S>'].mean()
        scale_factor = avg_sizeS
            
    #algo spec (for output)
    if(algo == 'EO' or algo == 'tauEO'):
        ALGOspec = (f"ND")
    if(algo == 'SA'):
        ALGOspec = (f"{param_name} = {param_value}, it = 2^{int( math.log(maxIters, 2.0) )}")
    if(algo == 'PA'):
        ALGOspec = (f"{param_name} = {param_value}, R = {replicas}, it = 2^{int(math.log(maxIters, 2.0))}")
    if(algo == 'CEA'):
        ALGOspec = (f"ND")
    if(algo == 'BP'):
        ALGOspec = (f"f = {f}, {param_name} = {param_value}")
    if(algo == 'BPR'):
        ALGOspec = (f"f = {f}, {param_name} = {param_value}, R = {replicas}")
    if(algo == 'BPA'):
        ALGOspec = (f"f = {f}, {param_name} = {param_value}")
    
    data[filename]['algo'] = algo + ': ' + ALGOspec
        
    #Set a musk according to the tolerance threshold for the results
    Results = Results.copy()
    mask = np.abs( ( round(Results['bestE'],4) - round(Egs,4) )/Egs ) <= tolerance
    Results.loc[mask, 'bestE'] = Egs
    data[filename]['dEmin'] = np.abs( ( round(Results['bestE'],4) - round(Egs,4) )/Egs )
    
    num_minE = Results[mask]['bestE'].count()
    success = num_minE / runs
    itersMinE = Results[ mask ]['bestIter'].mean()

    data[filename]['MinIter'] = Results['bestIter']

    return {
        'algo': algo,
        'algoSPEC': ALGOspec,
        'algoLABEL': algo + ': ' + ALGOspec,
        'param_name': param_name,
        'param_value': param_value,
        'f': f,
        '<f>': avg_f,
        '<BPit>': avg_BPit,
        'replicas': replicas,
        'bestE': bestE,
        '<|S|>': avg_sizeS,
        '<q>': avg_q,
        '<bestE>': avg_bestE,
        'err(bestE)': std_bestE,
        '<totIters>': avg_totIters,
        '<bestIters>': avg_bestIter,
        'scale_factor': scale_factor,
        'accepted': accepted,
        'success': success,
        'itersMinE': itersMinE,
        'runs': runs,
        'time': time
#        'spins': spins
    }


#Prepare data for analysis
def get_trajectories(dataType, data, INFO):

    if dataType == 'BP':
        trajectories = data
        trajectories['ID'] = (trajectories['iter'].diff() < 0).cumsum() + 1  # numbering runs from 1 to maxRuns
        trajectories = [group for _, group in trajectories.groupby('ID')]
        
    else:
        trajectories = []
        n = 1
        for file in data['filename'].unique():
            runs = INFO[file]['runs']
            for i in np.arange(0, runs):
                trajectory = data[data['filename'] == file].iloc[i::runs].copy()
                trajectory['ID'] = (i + 1) * n
                trajectory.rename(columns={"bestIter": "iter", "bestE": "E"}, inplace=True)
                trajectories.append(trajectory)
                i += 1
            n += 1
    
    return trajectories


#Compute dE vs T trajectory with bootstrap resampling
def analyze_trajectories(trajectories, B):

    Tmin = 1
    EbestEVOboot = []

    for _ in range(B):
        random.shuffle(trajectories)
        EVOshuffled = pd.concat(trajectories, ignore_index=True)
        EVOshuffled['iterTOT'] = EVOshuffled['iter'] * np.sort(EVOshuffled['ID'])
        EVOshuffled = EVOshuffled.sort_values('iterTOT').reset_index(drop=True)
        cummin_changes = EVOshuffled['dE'].cummin() != EVOshuffled['dE'].cummin().shift(1)
        cummin_EVO = EVOshuffled[cummin_changes]
        EbestEVOboot.append( cummin_EVO[['iterTOT', 'dE']] )

    # Expand to common time axis
    Tmax = max(df['iterTOT'].max() for df in EbestEVOboot)
    T1 = np.arange(64, 1000, 4)          # fine resolution at early times
    T2 = np.unique(np.round(np.logspace(np.log10(1000), np.log10(Tmax), 1000)).astype(int))
    Tfull = pd.DataFrame({'iterTOT': np.unique(np.concatenate([T1, T2]))})

    EbestEVOexp = []
    for df in EbestEVOboot:
        df_expanded = pd.merge_asof(Tfull, df, on='iterTOT', direction='backward')
        df_expanded['dE'] = df_expanded['dE'].fillna(method='ffill')
        df_expanded['dE'] = df_expanded['dE'].fillna(df['dE'].iloc[0])
        EbestEVOexp.append(df_expanded)
        
    # Compute average and quantiles
    dE_matrix = pd.concat([df['dE'] for df in EbestEVOexp], axis=1)
    dE_mean = dE_matrix.mean(axis=1)
    dE_std  = dE_matrix.std(axis=1)

    df_avg = pd.DataFrame({
        'T':        Tfull['iterTOT'].values,
        'dEavrg':   dE_mean,
        'dE_lower': np.maximum(dE_mean - dE_std, dE_matrix.min(axis=1)), #If dE_lower < 0. (which is nonsense), dE_lower = dE_min
        'dE_upper': dE_mean + dE_std,
    })
                
    return df_avg, EbestEVOboot


# Print info on file
def print_info(N, M, Egs, INFO, output):

    output.write( "#*************** ALGORITHMS STATS ***************\n\n" )
    output.write( "# Graph Info: N = {:d}, M = {:}\n".format(N, M) )
    for path, energy in Egs.items():
        output.write("# Best energy for {:>25s} Eg = {:.6f}\n".format(path, energy))
    output.write("\n# Samples Stats\n")
    output.write("#{:>5s}\t{:>5s}\t{:>35s}\t{:>10s}\t{:>8s}\t{:>13s}\t{:>12s}\t{:>12s}\t{:>12s}\t{:>6s}\t{:>5s}\t{:>5s}\t{:>5s}\t{:>7s}\t{:>5s}\t{:>9s}\n".format("Graph", "ALGO", "Parameters", "<Accepted>", "<bestIt>", "<bestIt(res)>", "Emin", "<Ebest>", "std(<Ebest>)", "<f>", "<|S|>", "<BPit>", "<q>", "Success", "Runs", "Time(s)"))

    for path, params in INFO.items():
        graphN = os.path.basename( os.path.dirname(path) ) if path else "ND"
        output.write("{:>6s}\t{:>5s}\t{:>35s}\t{:>10.3f}\t{:>8.2e}\t{:>13.2e}\t{:>12.6f}\t{:>12.3f}\t{:>12.3f}\t{:>6.1f}\t{:>5.1f}\t{:>6.1f}\t{:>5.1f}\t{:>7.3f}\t{:>5d}\t{:>9.1f}\n".format(graphN, params['algo'], params['algoSPEC'], params['accepted'] / params['<totIters>'], params['<bestIters>'],params['<bestIters>'] * params['scale_factor'], params['bestE'], params['<bestE>'], params['err(bestE)'],params['<f>'], params['<|S|>'], params['<BPit>'],params['<q>'],params['success'],params['runs'],params['time']))


# Plot algo performance
def plot_info(N, M, data, infoALL, output):

    infoALL = pd.DataFrame.from_dict(INFO, orient='index').reset_index().rename(columns={'index':'filename'})
    algorithms = sorted(infoALL['algoLABEL'].unique(), key=sort_ALGOlabels)
    infoALL['algoLABEL'] = pd.Categorical(infoALL['algoLABEL'], categories=algorithms, ordered=True)
    infoALL = infoALL.sort_values('algoLABEL')

    data = {k: v for k, v in data.items() if 'EVO' not in k}
    deltaEmin = pd.concat([pd.DataFrame({"algo": v['algo'], "dE": v['dEmin'], "iters": v['MinIter']})for k, v in data.items()], ignore_index=True).dropna()
    
    # Precompute metrics for top panel
    successAVRG = infoALL.groupby('algoLABEL', sort=False)['success'].mean()
    infoALL['<bestIters>scaled'] = infoALL['<bestIters>'] * infoALL['scale_factor']
    bestIt = infoALL.groupby('algoLABEL', sort=False)['<bestIters>scaled']
    bestItAVRG = infoALL.groupby('algoLABEL', sort=False)['<bestIters>scaled'].mean()
    bestItSTD = infoALL.groupby('algoLABEL', sort=False)['<bestIters>scaled'].std().fillna(0)

    output.write("\n# Averge Stats\n")
    output.write("#{:>36s}\t{:>10s}\t{:>14s}\t{:>12s}\t{:>13s}\t{:>7s}\n".format("ALGO", "<dEmin/Eb>", "Err(<dEmin/Eb>)", "<minIt>", "Err(<minIt>)", "Success"))
    deltaE_grouped = deltaEmin.groupby('algo')['dE']
    for algo in algorithms:
        output.write("{:>37s}\t{:>10.6f}\t{:>15.6f}\t{:>12.2e}\t{:>13.2e}\t{:>7.3f}\n".format(algo,deltaE_grouped.get_group(algo).mean(),deltaE_grouped.get_group(algo).std(),bestItAVRG.loc[algo],bestItSTD.loc[algo],successAVRG.loc[algo]))
    
    # Create figure with two panels
    fig, axs = plt.subplots( nrows=2, ncols=1, figsize=(7, 6.5), gridspec_kw={'height_ratios': [1, 2]}, sharex=True )
    fig.suptitle(f'Results Summary: N={N}, M={M}', fontsize=14)

    # --- Top panel: scatter of accepted_ratio ---
    ax1 = axs[0]
    ax1.grid(True, alpha=0.3)
    ax1.scatter( successAVRG.index, successAVRG.values, marker='o', s=120, ec='tab:blue', fc='tab:blue', lw=3.5, alpha=1, zorder=0)
    ax1.margins(x=0.15, y=0.2)
    ax1.set_ylabel('Success', color='tab:blue', fontsize=14)
    ax1.tick_params(axis='y', labelcolor='tab:blue')

    ax2 = ax1.twinx()
    ax2.errorbar(bestItAVRG.index, bestItAVRG, yerr=bestItSTD.replace([np.inf, -np.inf], 0), fmt='none', ecolor='tab:red', elinewidth=2, capsize=5, capthick=2, zorder=0)
    ax2.scatter( bestItAVRG.index, bestItAVRG, marker='^', s=100, ec='tab:red', fc='tab:red', lw=3, alpha=1, zorder=1)
    ax2.margins(x=0.15, y=0.2)
    ax2.set_ylabel('<Iters>', color='tab:red', fontsize=14)
    ax2.tick_params(axis='y', labelcolor='tab:red')
#    ax2.set_yscale('log')

    # --- Bottom panel: violin plot ---
    sns.violinplot(data=deltaEmin, x='algo', y='dE', order=algorithms, density_norm='width', common_norm=True, bw_adjust=.5, inner='box', alpha=0.7, cut=0, ax=axs[1])
    axs[1].axhline(0, color="black", linestyle="--", linewidth=1)
    axs[1].set_ylabel(r'$\Delta E_{min}  = (E_{best} - E_{min}) / E_{best}$', fontsize=14)
    axs[1].set_xlabel('Algorithm', fontsize=14)
    axs[1].grid(True, alpha=0.3)

    # Rotate tick labels properly
    axs[1].set_xticks(range(len(algorithms)))
    axs[1].set_xticklabels(algorithms, rotation=45, ha='right')

    plt.tight_layout(rect=[0, 0, 1, 0.96])  # leave space for suptitle
    plt.savefig(f'N{N}_M{M}_tol{tolerance}_results.pdf', format="pdf")
    plt.close(fig)


#Plot dE vs q (+ dE histo & q histo)
def plot_correlation(data, INFO):

    for file in data['filename'].apply(os.path.basename).unique():

        file_data = data[ data['filename'].apply(os.path.basename) == file ][['iter', 'itBP', 'dEi', 'dEb', 'q']]
        
        df_plot = file_data[ np.abs(file_data['dEi'] + file_data['dEb']) < 20.  ]
        dE = df_plot['dEi'] + df_plot['dEb']
        flips = df_plot['q']
        bp_iters = df_plot['itBP']
        dt = df_plot['iter'].diff().dropna()

        # --- Figure layout
        fig = plt.figure(figsize=(8, 8))
        gs = GridSpec(6, 6, figure=fig, wspace=0.05, hspace=0.05)

        # Subplots
        ax_yhist = fig.add_subplot(gs[1:6, 0])         # Left marginal
        ax_main = fig.add_subplot(gs[1:6, 1:6])        # Main scatter
        ax_xhist = fig.add_subplot(gs[0, 1:6], sharex=ax_main)  # Top marginal

        # --- Top histogram (flips)
        density = gaussian_kde(flips)
        x_range = np.linspace(flips.min(), flips.max(), 100)
        kde_values = density(x_range)
        kde_counts = kde_values * len(flips) * (flips.max() - flips.min()) / len(x_range)

        ax_xhist.plot(x_range, kde_counts, color='tab:blue', lw=4., alpha=1.)
        ax_xhist.fill_between(x_range, 0, kde_counts, alpha=0.3)
        ax_xhist.set_ylim(0.09, kde_counts.max() * 10)
        ax_xhist.set_yscale('log')

        ax_xhist.spines['top'].set_visible(True)
        ax_xhist.spines['right'].set_visible(True)
        ax_xhist.spines['left'].set_visible(True)
        ax_xhist.spines['bottom'].set_visible(True)
        ax_xhist.tick_params(bottom=False, labelbottom=False)
        ax_xhist.tick_params(axis='x', direction='in')  # Flip x ticks inward
        ax_xhist.tick_params(axis='y', which ='both', direction='in')  # Flip y ticks inward

        # --- Left histogram (dE)
        density = gaussian_kde(dE)
        y_range = np.linspace(dE.min(), dE.max(), 100)
        kde_values = density(y_range)
        kde_counts = kde_values * len(dE) * (dE.max() - dE.min()) / len(y_range)

        ax_yhist.plot(kde_counts, y_range, color='tab:blue', lw=4., alpha=1.)
        ax_yhist.fill_betweenx(y_range, 0, kde_counts, alpha=0.3)
        ax_yhist.margins(x=0.15, y=0.05)
        ax_yhist.set_xscale('log')
        
        ax_yhist.invert_xaxis()
        ax_yhist.spines['top'].set_visible(True)
        ax_yhist.spines['right'].set_visible(True)
        ax_yhist.spines['left'].set_visible(True)
        ax_yhist.spines['bottom'].set_visible(True)
        ax_yhist.tick_params(which ='both', direction = 'in', left=True, labelleft=True)
        ax_yhist.set_ylabel(r'$\Delta E_c = E_c(t) - E_c(t-1)$', labelpad=10)

        # --- Main scatter plot
        sc = ax_main.scatter(flips, dE, c=bp_iters, cmap='rainbow', vmin=0, vmax=100, s=40, edgecolor='white', linewidth=0.2, rasterized=True)
        ax_main.set_xlabel('Flipped Spins')
        ax_main.set_yticks([])
        ax_main.set_ylabel('')  # Moved to left marginal
        ax_main.tick_params(axis='x', direction='in')  # Flip x-ticks inward
        ax_main.tick_params(axis='y', direction='in')
        ax_main.margins(x=0.05, y=0.05)


        # --- Colorbar (vertical, outside grid, height-aligned)
        pos = ax_main.get_position()
        cb_ax = fig.add_axes([pos.x1 + 0.02, pos.y0, 0.02, pos.height])  # [left, bottom, width, height]
        cbar = fig.colorbar(sc, cax=cb_ax, orientation='vertical')
        cbar.set_label('BP iters')

        idx = [k for k in INFO.keys() if os.path.basename(k.replace('EVO', '')) == file.replace('EVO', '')]
        fig.suptitle(f"N = {N}, M = {M}\n{INFO[idx[0]]['algoLABEL']}", y=0.96)
        TitlefigCorr = Figtitle + f"f{INFO[idx[0]]['f']}_{INFO[idx[0]]['param_name']}_{INFO[idx[0]]['param_value']}_{INFO[idx[0]]['algo']}_Corr.pdf"
        plt.savefig(TitlefigCorr, format="pdf", bbox_inches='tight', dpi=300)


#Plot BP stats (cluster size, BP iters, acceptance ratio VS. ALGO)
def plot_BPstats(data, INFO):
    # Separate data sources
    base_files = data['filename'].apply(os.path.basename)
    data_noEVO = data[~base_files.str.contains('EVO')]
    data_EVO = data[base_files.str.contains('EVO')]
    
    # Prepare cluster size data
    conv_size = []
    for file in data_noEVO['filename'].apply(os.path.basename).unique():
        file_data = data_noEVO[data_noEVO['filename'].apply(os.path.basename) == file]
        maxIters = file_data['Iter'].iloc[-1]
        Results = file_data[file_data['Iter'] == maxIters]
        full_filename = file_data['filename'].iloc[0]
        f = INFO[full_filename]['<f>']
        for _, row in Results.iterrows():
            conv_size.append({'<|S|>': row['<|S|>'], 'algoLABEL': INFO[full_filename]['algoLABEL'], 'f': f})
    df_size = pd.DataFrame(conv_size)
    
    # Prepare acceptance rate data
    conv_accept = []
    for file in data_noEVO['filename'].apply(os.path.basename).unique():
        file_data = data_noEVO[data_noEVO['filename'].apply(os.path.basename) == file]
        maxIters = file_data['Iter'].iloc[-1]
        Results = file_data[file_data['Iter'] == maxIters]
        full_filename = file_data['filename'].iloc[0]
        f = INFO[full_filename]['<f>']
        for _, row in Results.iterrows():
            conv_accept.append({'accepted': row['accepted'] / maxIters, 'algoLABEL': INFO[full_filename]['algoLABEL'], 'f': f})
    df_accept = pd.DataFrame(conv_accept)
    
    # Prepare BP steps data
    data_list_bp = []
    for file in data_EVO['filename'].unique():
        filename = file.replace('EVO.bp', '.bp')
        algo_label = INFO[filename]['algoLABEL']
        bp_steps = data_EVO.loc[data_EVO['filename'] == file, 'itBP'].values
        data_list_bp.extend([{'<BPsteps>': bp, 'algoLABEL': algo_label, 'f': INFO[filename]['<f>']} for bp in bp_steps])
    df_bp = pd.DataFrame(data_list_bp)
    
    # Extract families and create common ordering
    for df in [df_size, df_accept, df_bp]:
        df['family'] = df['algoLABEL'].str.replace(r'\s+f\s*=\s*[\d.]+', '', regex=True).str.strip()
    
    all_labels = set(df_size['algoLABEL'].unique()) | set(df_accept['algoLABEL'].unique()) | set(df_bp['algoLABEL'].unique())
    label_order = sorted(all_labels, key=lambda s: [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', s)])
    
    all_families = set(df_size['family'].unique()) | set(df_accept['family'].unique()) | set(df_bp['family'].unique())
    family_colors = dict(zip(all_families, sns.color_palette('Set2', n_colors=len(all_families))))
    single_family = len(all_families) == 1
    
    # Define positions and common properties
    label_positions = {l: df_size[df_size['algoLABEL']==l]['f'].iloc[0] if single_family else label_order.index(l) for l in label_order}
    common_colors = [family_colors[df_size[df_size['algoLABEL']==l]['family'].iloc[0]] for l in label_order if l in df_size['algoLABEL'].values]
    common_labels = [l for l in label_order if l in df_size['algoLABEL'].values]
    
    # Helper function for boxplot
    def draw_boxplot(ax, df, y_col, data_labels):
        if single_family:
            for label in data_labels:
                f_val = label_positions[label]
                values = df[df['algoLABEL'] == label][y_col].values
                ax.boxplot([values], positions=[f_val], widths=[f_val * 0.3], patch_artist=True, showfliers=False,
                          capwidths=[f_val * 0.2], capprops=dict(color='dimgrey', lw=1.5),
                          whiskerprops=dict(color='dimgrey', lw=1.5),
                          boxprops=dict(color='dimgrey', facecolor=common_colors[data_labels.index(label)], lw=1.5),
                          medianprops=dict(color='dimgrey', lw=1.5))
            if len(data_labels) > 1:
                fam_labels = sorted(data_labels, key=lambda l: label_positions[l])
                pos = [label_positions[l] for l in fam_labels]
                med = [df[df['algoLABEL']==l][y_col].median() for l in fam_labels]
                ax.plot(pos, med, 'o-', color=family_colors[list(df['family'].unique())[0]], linewidth=3, markersize=0.1, alpha=1., zorder=0)
        else:
            sns.boxplot(x='algoLABEL', y=y_col, data=df, width=0.65, showfliers=False, hue='algoLABEL',
                       palette=common_colors, order=data_labels, ax=ax, legend=False)
            for family in df['family'].unique():
                fam_labels = [l for l in data_labels if df[df['algoLABEL']==l]['family'].iloc[0] == family]
                if len(fam_labels) > 1:
                    pos = [label_positions[l] for l in fam_labels]
                    med = [df[df['algoLABEL']==l][y_col].median() for l in fam_labels]
                    ax.plot(pos, med, 'o-', color=family_colors[family], linewidth=3, markersize=0.1, alpha=1., zorder=0)
    
    # Plot
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    
    # Subplot 1: Cluster Size
    for family in df_size['family'].unique():
        fam_labels = [l for l in label_order if l in df_size['algoLABEL'].values and df_size[df_size['algoLABEL']==l]['family'].iloc[0] == family]
        if fam_labels:
            pos = [label_positions[l] for l in fam_labels]
            means = [df_size[df_size['algoLABEL']==l]['<|S|>'].mean() for l in fam_labels]
            stds = [df_size[df_size['algoLABEL']==l]['<|S|>'].std() for l in fam_labels]
            ax1.errorbar(pos, means, yerr=stds, fmt='o-' if single_family else 'o', color=family_colors[family],
                        linewidth=3, markersize=10, capsize=5, capthick=2, alpha=1.)
    ax1.set_ylabel("Cluster Size", fontsize=14)
#    ax1.set_ylim(10, 7000)
    ax1.set_yscale('log')
    
    # Subplot 2: Acceptance Rate
    draw_boxplot(ax2, df_accept, 'accepted', common_labels)
    ax2.set_ylabel("Acceptance Rate", fontsize=14)
    ax2.yaxis.set_major_locator(MultipleLocator(0.2))
    ax2.set_ylim(-0.05, 1.05)
    ax2.set_xlabel('')
    
    # Subplot 3: BP Steps
    draw_boxplot(ax3, df_bp, '<BPsteps>', common_labels)
    ax3.set_ylabel("BP Iters", fontsize=14)
    ax3.tick_params(direction='in', length=4, width=1.5)
    ax3.yaxis.set_major_locator(MultipleLocator(20))
    ax3.set_ylim(0.7, 150)
    ax3.set_yscale('log')
    
    # Set x-axis
    if single_family:
        ax3.set_xlabel("Frustration", fontsize=14)
#        ax3.set_xlim(0.75, 2000)
        ax3.set_xscale('log')
    else:
        ax3.set_xlabel("Algorithm", fontsize=14)
        ax3.set_xticks(range(len(label_order)))
        ax3.set_xticklabels(label_order, rotation=45, ha='right')
    
    plt.suptitle(f'N = {N}, M = {M}', fontsize=14)
    plt.tight_layout()
    plt.savefig(Figtitle + '_BPstats.pdf', format="pdf")
    
    
    
#Plot histogram of dIters
def plot_iterations(data, INFO):

    plt.figure(figsize=(8, 6))
    i = 0
    for file in data['filename'].apply(os.path.basename).unique():
        file_data = data[ data['filename'].apply(os.path.basename) == file ]
        dt = file_data['iter'].diff().dropna()
        dt = dt[ dt > 0. ]
        full_filename = file_data['filename'].iloc[0].replace('EVO', '')
        plt.hist(dt, bins='rice', log=True, histtype='step', ec=clrBP[i], lw=3., alpha=0.9, label=INFO[full_filename]['algoLABEL'])
        i += 1

    plt.tick_params(direction='in', length=4, width=1.5)
    plt.xticks(rotation=45, ha='right')
    plt.title('N = ' + str(N) + ', M = ' + M, fontsize=14)
    plt.xlabel("$\Delta t$", fontsize=14)
    plt.ylabel("Count", fontsize=14)
    plt.xscale('log')
    plt.legend()
    plt.tight_layout(pad=3.0)
    plt.savefig(Figtitle + '_BPiters.pdf', format="pdf")
    


# ----------------------------------------------------------
# ------------------------ MAIN ----------------------------

if len(sys.argv) < 4:
    sys.stderr.write("\nUsage: {} OUTPUT tolerance <data files>\n".format(sys.argv[0]))
    sys.stderr.write("OUTPUT options: FULL, OVERVIEW, INFO\n")
    sys.stderr.write("\t> FULL: All figures and data files\n")
    sys.stderr.write("\t> OVERVIEW: Only EbestEVO figure and data files\n")
    sys.stderr.write("\t> INFO: Only data file\n\n")
    sys.exit(1)

OUTPUT = sys.argv[1].upper()
if OUTPUT not in ['FULL', 'OVERVIEW', 'INFO']:
    sys.stderr.write("Error: OUTPUT must be FULL, OVERVIEW, or INFO\n")
    sys.exit(1)

tolerance = float(sys.argv[2])
infiles = sys.argv[3:]

#Get N and M from file name
first_filename = os.path.basename(infiles[0])
N = int( re.search('N(.*?)_', first_filename).group(1) )
M = re.search('_M(.*?)_', first_filename).group(1)

#Init plots
Figtitle = 'N' + str(N) + '_M' + M + '_tol' + str(tolerance)

figEbestEVO = None
if OUTPUT in ['FULL', 'OVERVIEW']:
    figEbestEVO = plt.figure(figsize=(7, 5)) #Use (7, 5) for paper figs, (7, 10) for slides figs.

# Get best energies: group files by directory and find min energy for each group
Egs = {}
infiles_sorted = sorted(infiles, key=os.path.dirname)
for path, files_group in groupby(infiles_sorted, key=os.path.dirname):
    files_in_path = list(files_group)
    energies = extract_optimal_energies(files_in_path)
    Egs[path] = np.min(energies)

#Import all data
data = {file: import_data(file) for file in infiles}

INFO ={}
for file in infiles:
    if 'EVO' not in file:
        INFO[file] = get_info(file, data, Egs[os.path.dirname(file)], N)

#Order data according to algo (type & specs)
def sort_ALGOlabels(label):
    prefix = next((i for i, s in enumerate(['EO', 'tauEO', 'SA', 'PA', 'CEA']) if s in label.upper()), 999)
    natural_key = [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', label)]
    return (prefix == 999, prefix, natural_key)
    
Egs = dict(sorted(Egs.items(), key=lambda x: sort_ALGOlabels(x[0])))
data = dict(sorted(data.items(), key=lambda x: sort_ALGOlabels(x[0])))
INFO = dict(sorted(INFO.items(), key=lambda x: sort_ALGOlabels(x[0])))

#Rename dictionaries (data, INFO) keys to avoid filename discrepancies
data = {
    (os.path.join(os.path.dirname(key), re.sub(r"N[^_]*_M[^_]*_", "", os.path.basename(key)))
     if os.path.dirname(key) else re.sub(r"N[^_]*_M[^_]*_", "", os.path.basename(key))): value
    for key, value in data.items()
}

INFO = {
    (os.path.join(os.path.dirname(key), re.sub(r"N[^_]*_M[^_]*_", "", os.path.basename(key)))
     if os.path.dirname(key) else re.sub(r"N[^_]*_M[^_]*_", "", os.path.basename(key))): value
    for key, value in INFO.items()
}


#Default output (INFO)
output = open('N' + str(N) + '_M' + M + '_tol' + str(tolerance) + '_results.dat', 'w+')
print_info(N, M, Egs, INFO, output)
plot_info(N, M, data, INFO, output)
output.close()


# ***** Data Analisys ****

# ***** EO/tauEO *****
eo_variants = [
    ('.eo',    'EO',    'green'),
    ('.taueo', 'tauEO', 'lightgreen'),
]

for ext, algo_label, color in eo_variants:
    if not any(ext in filename for filename in data.keys()):
        continue

    eodata = {k: v for k, v in data.items() if ext in k}
    for filename in eodata.keys():
        Egs_val = Egs[os.path.dirname(filename)]
        mask = np.fabs((round(eodata[filename]['bestE'], 4) - round(Egs_val, 4)) / Egs_val) <= tolerance
        eodata[filename].loc[mask, 'bestE'] = Egs_val
        eodata[filename]['dE'] = eodata[filename]['bestE'] - Egs_val

    dataALL = pd.concat(eodata.values(), keys=eodata.keys()).reset_index(level=0).rename(columns={'level_0': 'filename'})

    if OUTPUT in ['OVERVIEW', 'FULL']:
        for file in dataALL['filename'].apply(os.path.basename).unique():
            trajectories = get_trajectories(algo_label, dataALL[dataALL['filename'].apply(os.path.basename) == file][['filename', 'bestIter', 'bestE', 'dE']], INFO)
            df_avg, EbestEVOexp = analyze_trajectories(trajectories, B)

            matching_keys = [k for k in INFO.keys() if os.path.basename(k) == file]
            label = f"{INFO[matching_keys[0]]['algoLABEL']}"
            plt.figure(figEbestEVO.number)
            plt.plot(df_avg['T'], df_avg['dEavrg'], color=color, label=label)
            plt.fill_between(df_avg['T'], df_avg['dE_lower'], df_avg['dE_upper'], color=color, alpha=0.3)

# ***** SA *****
if any('.sa' in filename for filename in data.keys()):

    i = 0
    SAdata = {k: v for k, v in data.items() if '.sa' in k}
    for filename in SAdata.keys():
        Egs_val = Egs[os.path.dirname(filename)]
        mask = np.fabs((round(SAdata[filename]['bestE'],4) - round(Egs_val,4)) / Egs_val) <= tolerance
        SAdata[filename].loc[mask, 'bestE'] = Egs_val
        SAdata[filename]['dE'] = round(SAdata[filename]['bestE'],4) - round(Egs_val,4)
        SAdata[filename]['dE'] = SAdata[filename]['bestE'] - Egs_val
        
    dataALL = pd.concat(SAdata.values(), keys=SAdata.keys()).reset_index(level=0).rename(columns={'level_0': 'filename'})

    if OUTPUT in ['OVERVIEW', 'FULL']:

        for file in dataALL['filename'].apply(os.path.basename).unique():
            trajectories = get_trajectories('SA', dataALL[ dataALL['filename'].apply(os.path.basename) == file ][['filename', 'bestIter','bestE','dE']], INFO )
            df_avg, EbestEVOexp = analyze_trajectories(trajectories, B)

            matching_keys = [k for k in INFO.keys() if os.path.basename(k) == file]
            label = f"{INFO[matching_keys[0]]['algoLABEL']}"
            plt.figure(figEbestEVO.number)
            plt.plot( N * df_avg['T'], df_avg['dEavrg'], color = clrMC[i], label = label ) #Avrg scale factor = N
            plt.fill_between( N * df_avg['T'], df_avg['dE_lower'], df_avg['dE_upper'], color = clrMC[i], alpha=0.3)
            
            i += 1
            
# ***** PA *****
if any('.pa' in filename for filename in data.keys()):

    i = 0
    PAdata = {k: v for k, v in data.items() if '.pa' in k}
    for filename in PAdata.keys():
        Egs_val = Egs[os.path.dirname(filename)]
        mask = np.fabs((round(PAdata[filename]['bestE'],4) - round(Egs_val,4)) / Egs_val) <= tolerance
        PAdata[filename].loc[mask, 'bestE'] = Egs_val
        PAdata[filename]['dE'] = round(PAdata[filename]['bestE'],4) - round(Egs_val,4)
        PAdata[filename]['dE'] = PAdata[filename]['bestE'] - Egs_val
        
    dataALL = pd.concat(PAdata.values(), keys=PAdata.keys()).reset_index(level=0).rename(columns={'level_0': 'filename'})

    if OUTPUT in ['OVERVIEW', 'FULL']:

        for file in dataALL['filename'].apply(os.path.basename).unique():
            trajectories = get_trajectories('PA', dataALL[ dataALL['filename'].apply(os.path.basename) == file ][['filename', 'bestIter','bestE','dE']], INFO )
            df_avg, EbestEVOexp = analyze_trajectories(trajectories, B)

            matching_keys = [k for k in INFO.keys() if os.path.basename(k) == file]
            label = f"{INFO[matching_keys[0]]['algoLABEL']}"
            plt.figure(figEbestEVO.number)
            plt.plot( INFO[matching_keys[0]]['scale_factor'] * df_avg['T'], df_avg['dEavrg'], color = clrMC[i+3], label = label ) #Avrg scale factor = N
            plt.fill_between( INFO[matching_keys[0]]['scale_factor'] * df_avg['T'], df_avg['dE_lower'], df_avg['dE_upper'], color = clrMC[i+3], alpha=0.3)
            
            i += 1

# ***** CEA *****
if any('.cea' in filename for filename in data.keys()):

    CEAdata = {k: v for k, v in data.items() if '.cea' in k}
    for filename in CEAdata.keys():
        Egs_val = Egs[os.path.dirname(filename)]
        CEAdata[filename] = CEAdata[filename].copy()
        mask = np.fabs((round(CEAdata[filename]['bestE'],4) - round(Egs_val,4)) / Egs_val) <= tolerance
        CEAdata[filename].loc[mask, 'bestE'] = Egs_val
        CEAdata[filename]['dE'] = round(CEAdata[filename]['bestE'],4) - round(Egs_val,4)
    dataALL = pd.concat(CEAdata.values(), keys=CEAdata.keys()).reset_index(level=0).rename(columns={'level_0': 'filename'})

    if OUTPUT in ['OVERVIEW', 'FULL']:

        for file in dataALL['filename'].apply(os.path.basename).unique():
            trajectories = get_trajectories('CEA', dataALL[ dataALL['filename'].apply(os.path.basename) == file ][['filename','bestIter','bestE','dE']], INFO )
            df_avg, EbestEVOexp = analyze_trajectories(trajectories, B)

            #Plot on 'figEbestEVO'
            matching_keys = [k for k in INFO.keys() if os.path.basename(k) == file]
            
            #Compute avrg scale factor
            scale_factor_avg = 0.
            for k in matching_keys:
                scale_factor_avg += INFO[k]['scale_factor']
            scale_factor_avg /= len(matching_keys)
                
            label = f"{INFO[matching_keys[0]]['algoLABEL']}"
            plt.figure(figEbestEVO.number)
            plt.plot( scale_factor_avg * df_avg['T'], df_avg['dEavrg'], color = 'orange', label = label )
            plt.fill_between( df_avg['T'], df_avg['dE_lower'], df_avg['dE_upper'], color = 'orange', alpha=0.3)


# ***** BP *****
if any('.bp' in filename for filename in data.keys()):

    BPdata = {k: v for k, v in data.items() if '.bp' in k}

    dataALL = pd.concat(BPdata.values(), keys=BPdata.keys()).reset_index(level=0).rename(columns={'level_0': 'filename'})
    plot_BPstats(dataALL, INFO)

    #Select 'evo' data files
    if any('EVO' in filename for filename in BPdata.keys()):

        i = 0
        evo_data = {k: v for k, v in BPdata.items() if 'EVO' in k}
        for filename in evo_data.keys():
            Egs_val = Egs[os.path.dirname(filename)]
            evo_data[filename] = evo_data[filename].copy()
            mask = np.fabs((round(evo_data[filename]['E'],4) - round(Egs_val,4)) / Egs_val) <= tolerance
            evo_data[filename].loc[mask, 'E'] = Egs_val
            evo_data[filename]['dE'] = round(evo_data[filename]['E'],4) - round(Egs_val,4)
        
        dataALL = pd.concat(evo_data.values(), keys=evo_data.keys()).reset_index(level=0).rename(columns={'level_0': 'filename'})

        if OUTPUT in ['OVERVIEW', 'FULL']:

            for file in dataALL['filename'].apply(os.path.basename).unique():
                # Note that the scale factor is only computed from the first .bp file as averages are more stable
                trajectories = get_trajectories('BP', dataALL[ dataALL['filename'].apply(os.path.basename) == file ][['iter','E','dE']].copy(), INFO)
                df_avg, EbestEVOexp = analyze_trajectories(trajectories, B)
                
                #Plot on 'figEbestEVO'
                matching_keys = [k for k in INFO.keys() if os.path.basename(k.replace('EVO', '')) == file.replace('EVO', '')]
                
                #Compute avrg scale factor
                scale_factor_avg = 0.
                for k in matching_keys:
                    scale_factor_avg += INFO[k]['scale_factor']
                scale_factor_avg /= len(matching_keys)
                                
                label = f"{INFO[matching_keys[0]]['algoLABEL']}"
                plt.figure(figEbestEVO.number)
                plt.plot( scale_factor_avg * df_avg['T'], df_avg['dEavrg'], color = clrBP[i], label = label )
                plt.fill_between( scale_factor_avg * df_avg['T'], df_avg['dE_lower'], df_avg['dE_upper'], color = clrBP[i], alpha=0.3)
                
                i += 1

        if OUTPUT == 'FULL':
            plot_iterations(dataALL, INFO)
            plot_correlation(dataALL, INFO)

    

# Print summary of what was generated
print(f"\nOutput mode: {OUTPUT}")
print("Generated:")
print("- Data file: N{}_M{}_tol{}_results.dat".format(N, M, tolerance))
print("- Figure: N{}_M{}_tol{}_results.pdf".format(N, M, tolerance))

if OUTPUT in ['OVERVIEW', 'FULL']:
    plt.figure(figEbestEVO.number)
    plt.xscale("log")
    plt.yscale("log")
    plt.title('N = ' + str(N) + ', M = ' + M)
    plt.xlabel('Elementary Ops. ($n$) per run')
    plt.ylabel('$\Delta E(n) = E_{best} - E(n)$')
    plt.legend()
    plt.savefig(Figtitle + '_EbestEVO.pdf', format="pdf")

    print("- Figure: {}_tol{}_EbestEVO.pdf".format(Figtitle[:-1], tolerance))

if OUTPUT == 'FULL':
        print("- Figure: {}_BPstats.pdf".format(Figtitle[:-1]))
        print("- Figure: {}_BPiters.pdf".format(Figtitle[:-1]))
        print("- Energy correlation figures for each BP EVO file")
print("\n")

sys.exit()
