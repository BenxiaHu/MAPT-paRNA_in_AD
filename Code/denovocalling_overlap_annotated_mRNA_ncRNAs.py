import matplotlib.pyplot as plt
labels = 'mRNA', 'lincRNA','antisenseRNA','other-ncRNAs','unknown-ncRNAs'
sizes = [12815, 1517,1191,1769,27127]

fig, ax = plt.subplots()

def autopct_format(values):
    def my_format(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{v:d}'.format(v=val)
    return my_format

ax.pie(sizes, labels = labels, autopct = autopct_format(sizes))
plt.savefig('denovocalling_overlap_annotated_mRNA_ncRNA_1.pdf')


labels = 'mRNA', 'lincRNA','antisenseRNA','other-ncRNAs','unknown-ncRNAs'
sizes = [12815, 1517,1191,1769,27127]

fig, ax = plt.subplots()
ax.pie(sizes, labels = labels,autopct='%1.1f%%')
plt.savefig('denovocalling_overlap_annotated_mRNA_ncRNA_2.pdf')
