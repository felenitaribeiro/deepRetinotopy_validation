import sys
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
import os
import numpy as np


# function to calculate Cohen's d for independent samples
def cohend(d1, d2):
	# calculate the size of samples
	n1, n2 = len(d1), len(d2)
	# calculate the variance of the samples
	s1, s2 = np.var(d1, ddof=1), np.var(d2, ddof=1)
	# calculate the pooled standard deviation
	s = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
	# calculate the means of the samples
	u1, u2 = np.mean(d1), np.mean(d2)
	# calculate the effect size
	return (u1 - u2) / s

def generate_boxplot_group_difference(method, results, output_dir):
    results_method = results[results['method'] == method].copy()
    
    # Create figure with 2x2 subplots
    fig, axes = plt.subplots(1, 2, figsize=(8, 5))

    # Plot 1: Cortical VMA (top-left)
    plt.sca(axes[1])
    plt.hlines(0, -.5, 1.5, color='black')
    if method == 'empirical':
        plt.ylim(-30, 150)
    else:
        plt.ylim(-20, 105)
    
    sns.barplot(data=results_method, y='cortical_vma', x="group", errorbar='se', ax=axes[1], palette=['#999999ff', '#ecececff'], linewidth= 3, width=0.4, hue='group', legend=False)
    sns.swarmplot(data=results_method, y='cortical_vma', x="group", ax=axes[1], color='black', size = 1)
    plt.ylabel('Asymmetry index', fontsize=17)
    plt.xlabel('')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=14)
    ttest_vma = scipy.stats.ttest_ind(results_method['cortical_vma'][results_method['group']=='adults'], results_method['cortical_vma'][results_method['group']=='children'])
    cohend_vma = cohend(results_method['cortical_vma'][results_method['group']=='adults'], results_method['cortical_vma'][results_method['group']=='children'])
    print(f'T-test for {method} VMA - p-value {ttest_vma.pvalue}, tstat {ttest_vma.statistic}')
    print(f"Cohen's d for {method} VMA: {cohend_vma}")
    print(f'Degree of freedom for HVA t-test: {len(results_method["cortical_hva"][results_method["group"]=="adults"]) + len(results_method["cortical_hva"][results_method["group"]=="children"]) - 2}')

    axes[1].set_title(f'Cortical VMA', fontsize=17)
    if ttest_vma.pvalue < 0.05:
        if ttest_vma.pvalue < 0.001:
            pvalue = '< 0.001'
            plt.text(0.5, 0.9, f'p {pvalue}', transform=axes[1].transAxes, ha='center', fontsize=14, color='red')

    # Plot 2: Cortical HVA (top-left)
    plt.sca(axes[0])
    plt.hlines(0, -.5, 1.5, color='black')
    if method == 'empirical':
        plt.ylim(-30, 150)
    else:
        plt.ylim(-20, 105)
    sns.barplot(data=results_method, y='cortical_hva', x="group", errorbar='se', ax=axes[0], palette=['#999999ff', '#ecececff'], linewidth= 3, width=0.4, hue='group', legend=False)
    sns.swarmplot(data=results_method, y='cortical_hva', x="group", ax=axes[0], color='black', size = 1)
    plt.ylabel('Asymmetry index', fontsize=17)
    plt.xlabel('')
    axes[0].set_title(f'Cortical HVA', fontsize=17)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=14)
    ttest_hva = scipy.stats.ttest_ind(results_method['cortical_hva'][results_method['group']=='adults'], results_method['cortical_hva'][results_method['group']=='children'])
    cohend_hva = cohend(results_method['cortical_hva'][results_method['group']=='adults'], results_method['cortical_hva'][results_method['group']=='children'])
    print(f'T-test for {method} HVA - p-value {ttest_hva.pvalue}, tstat {ttest_hva.statistic}')
    print(f"Cohen's d for {method} HVA: {cohend_hva}")
    print(f'Degree of freedom for HVA t-test: {len(results_method["cortical_hva"][results_method["group"]=="adults"]) + len(results_method["cortical_hva"][results_method["group"]=="children"]) - 2}')
    if ttest_hva.pvalue < 0.05:
        if ttest_hva.pvalue < 0.001:
            pvalue = '< 0.001'
            plt.text(0.5, 0.9, f'p {pvalue}', transform=axes[0].transAxes, ha='center', fontsize=14, color='red')

    if os.path.exists(output_dir) is False:
        os.makedirs(output_dir)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{output_dir}/barplot_group_difference_{method}_wedge45.pdf', dpi=300, bbox_inches='tight')

def main():
    parser = argparse.ArgumentParser(description='Visual Field Sampling Plotting Script')
    parser.add_argument('--results_abcd_csv', type=str, required=True, help='Path to the CSV file containing results data')
    parser.add_argument('--results_hcp_csv', type=str, required=True, help='Path to the CSV file containing HCP results data')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the output plots')
    args = parser.parse_args()


    results_abcd = pd.read_csv(args.results_abcd_csv)
    results_abcd['group'] = ['children'] * len(results_abcd)

    results_hcp = pd.read_csv(args.results_hcp_csv)
    results_hcp['group'] = ['adults'] * len(results_hcp)

    results = pd.concat([results_abcd, results_hcp], ignore_index=True)

    for method in ['deepRetinotopy']:
        generate_boxplot_group_difference(method, results, output_dir=args.output_dir)

if __name__ == '__main__':
    main()
