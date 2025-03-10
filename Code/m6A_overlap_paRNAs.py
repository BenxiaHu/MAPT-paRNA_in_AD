import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"
labels = 'm6A-paRNAs', 'non-m6A-paRNAs'
sizes = [1172, 1870]

fig, ax = plt.subplots()

def autopct_format(values):
    def my_format(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{v:d}'.format(v=val)
    return my_format

ax.pie(sizes, labels = labels, autopct = autopct_format(sizes))
plt.savefig('m6A-overlap-paRNAs.pdf')


labels = 'hyper-m6A-paRNAs','hypo-m6A-paRNAs', 'stable-m6A-paRNAs'
sizes = [91, 59,1022]

fig, ax = plt.subplots()

def autopct_format(values):
    def my_format(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{v:d}'.format(v=val)
    return my_format

ax.pie(sizes, labels = labels, autopct = autopct_format(sizes))
plt.savefig('m6A-overlap-paRNAs2.pdf')

