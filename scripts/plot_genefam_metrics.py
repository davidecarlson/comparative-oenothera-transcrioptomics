#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

oen_df = pd.read_csv('comp_trans_oen.csv', sep = ',', header = 0)

print(oen_df)

def r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2

def regressionplot(x, y, df, ci, xlim1, xlim2, xlab, ylab, png):
	slope, intercept, r_value, p_value, std_err = np.round(stats.linregress(df[x],df[y]), 3)
	sns.set(rc={'figure.figsize':(10,6)})
	sns.set_style("white")
	sns.set_context("paper", font_scale=1.25)
	sns.despine()
	ax = sns.regplot(x=x, y=y, data=df, ci=ci, fit_reg = True,
		line_kws={'label':f'r2 = {np.round(r_value**2, 3)}; p = {p_value}'})
		#line_kws={'label':"y={0:.1f}x+{1:.1f}".format(slope,intercept)})
	ax.legend()
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.xlim(xlim1, xlim2)
	plt.savefig(png, bbox_inches='tight')
	plt.show()
	plt.clf()

regressionplot('N50', 'Genes gained', oen_df, 95, 750, 1350, 'Transcript N50', 'No. of genes gained', 'n50_v_genegain.png')

regressionplot('No. Transcripts (longest isoform)', 'Genes gained', oen_df, 95, 32000, 50000, 'No. of Transcripts', 'No. of genes gained', 'numtrans_v_genegain.png')

regressionplot('Complete (%)', 'Genes gained', oen_df, 95, 40, 82, '% Complete BUSCOs', 'No. of genes gained', 'busco_v_genegain.png')

regressionplot('Complete (%)', 'Avg. Expansion', oen_df, 95, 40, 82, '% Complete BUSCOs', 'Avg. Expansion', 'busco_v_avgexpan.png')

regressionplot('Missing (%)', 'Genes lost', oen_df, 95, 12, 40, '% Missing BUSCOs', 'No. of genes lost', 'busco_v_geneloss.png')

