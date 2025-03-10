import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"
labels = 'down_MAPT_paRNA', 'other_down', 'up'
sizes = [62, 118, 23]

fig, ax = plt.subplots()

def autopct_format(values):
    def my_format(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{v:d}'.format(v=val)
    return my_format

ax.pie(sizes, labels = labels, autopct = autopct_format(sizes))
plt.savefig('MAPT_paRNA_DEGs_overlap_MUSCI.pdf')