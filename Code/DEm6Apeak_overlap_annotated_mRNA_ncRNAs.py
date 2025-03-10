import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np

def pieplot(sampleid):
    ncRNA = pd.read_csv("m6A_overlap_denovo_RNA2.txt",sep="\t",header=0)
    ncRNA = ncRNA[['type',sampleid]]
    labels = ncRNA['type']
    ncRNAsizes = ncRNA[sampleid]
    fig, ax = plt.subplots()
    def func(pct, allvals):
        absolute = int(np.round(pct/100.*np.sum(allvals)))
        return f"{pct:.1f}%\n({absolute:d})"
    wedges, texts, autotexts = ax.pie(ncRNAsizes, autopct=lambda pct: func(pct, ncRNAsizes),
                                  textprops=dict(color="w"))
    ax.legend(wedges, labels,
          title="type",
          loc="best",
          bbox_to_anchor=(1, 0, 0.5, 1))
    plt.setp(autotexts, size=8, weight="bold")
    ax.set_title(sampleid+'_annotated_unknown_ncRNAs_overlap_m6Apeak')
    plt.savefig(sampleid+'_m6Apeak_overlap_annotated_unknown_ncRNAs2.pdf')


pieplot("AD")
pieplot("Normal")
