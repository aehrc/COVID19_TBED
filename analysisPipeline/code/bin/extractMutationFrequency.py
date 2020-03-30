import sys
import numpy as np
import collections
import math
# import plotly.graph_objects as go
import psutil

#import data
data = {}
with open(sys.argv[1]) as f: 
    for line in f:
        line =  line.rstrip()
        if line[0] == ">":
            name = line
            if name not in data:
                data[name] = ""
        else:
            data[name] = data[name] + str(line)

list_data = [list(data[i]) for i in data]
array_data = np.array(list_data)

out1 = ''.join((sys.argv[2],'.uncondensedMutationFreq.txt'))
out2 = ''.join((sys.argv[2],'.condensedMutationFreq.txt'))

#make mutation frequency plot, includes Ns and deletions
values = []
num_strains = len(data)
for column in array_data.T: #find number of strains not matching consensus base
    values.append((num_strains - max(collections.Counter(column).values()))/float(num_strains))

with open(out1, 'w') as outfile1:
    for listitem in values:
        outfile1.write('%s\n' % listitem)

# #make x axis 
# random_x = np.linspace(0, len(values)-1, len(values))

# #make figure
# fig = go.Figure()
# fig.add_trace(go.Scatter(x=random_x, y=[math.sqrt(x) for x in values],
#                     mode='lines'))
# fig.update_xaxes(title_text="Genomic Coordinates (bp)")
# fig.update_yaxes(title_text="Mutation Frequency")
# fig.update_yaxes(range=[0, 0.6])

# fig.update_layout(
#     width=1600,
#     height=400,
#     title="Mutation Frequency across strains with Ns and deletions included")

# fig.update_layout(
#     xaxis = dict(
#         tickmode = 'linear',
#         tick0 = 0,
#         dtick = 3000)
#     )

# fig1 = ''.join((sys.argv[2],'IncludesDeletions.pdf'))
# fig2 = ''.join((sys.argv[2],'IncludesDeletionsSquished.pdf'))
# fig3 = ''.join((sys.argv[2],'ExcludesDeletions.pdf'))
# fig4 = ''.join((sys.argv[2],'ExcludesDeletionsSquished.pdf'))

# fig.write_image(fig1, width=1000, height=600)
# fig.write_image(fig2, width=1000, height=300)

#make mutation frequency plot, excludes Ns and deletions
values = []
num_strains = len(data)
for column in array_data.T:
    current = collections.Counter(column)
    if "N" in current:
        del current['N']
    if "-" in current:
        del current['-']
    maxCount = sum(current.values())
    values.append((maxCount - max(current.values()))/float(maxCount))

with open(out2, 'w') as outfile2:
    for listitem in values:
        outfile2.write('%s\n' % listitem)

# #make x axis
# random_x = np.linspace(0, len(values)-1, len(values))

# #make figure
# fig = go.Figure()
# fig.add_trace(go.Scatter(x=random_x, y=[math.sqrt(x) for x in values],
#                     mode='lines'))
# fig.update_xaxes(title_text="Genomic Coordinates (bp)")
# fig.update_yaxes(title_text="Mutation Frequency")
# fig.update_yaxes(range=[0, 0.6])
# fig.update_layout(
#     width=1600,
#     height=400,
#     title="Mutation Rate across strains with Ns and deletions removed")

# fig.update_layout(
#     xaxis = dict(
#         tickmode = 'linear',
#         tick0 = 0,
#         dtick = 3000)
#     )

# fig.write_image(fig3, width=1000, height=600)
# fig.write_image(fig4, width=1000, height=300)

