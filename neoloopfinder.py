from neoloop.visualize.core import *
import cooler
import pandas as pd
import matplotlib.pyplot as plt

sample = '236465' # HiC sample id

clr = cooler.Cooler('data_HOXC13/'+sample+'-HC01_ICE.mcool::resolutions/5000')
assembly = 'C8	deletion,12,20270000,+,12,53900000,-	12,19700000	12,54005000' # select assembly to plot
loops = 'data_HOXC13/'+sample+'-HC01.neo-loops.txt' # load neoloop data
outdir = 'plots'
List = [line.rstrip() for line in open('data/allOnco-genes.txt')] # load oncogenes list

vis = Triangle(clr, assembly, n_rows=6, figsize=(40, 36), track_partition=[5, 0.4, 0.8, 0.8, 0.8, 0.5], span=500000, correct='weight', space=0.03)
vis.matrix_plot(vmin=0, cbr_fontsize=40)
vis.plot_chromosome_bounds(linewidth=6)
vis.plot_loops(loops, face_color='none', edgecolors='#000000', marker_type='o', marker_size=600, cluster=True, onlyneo=True) # only show neo-loops
vis.plot_genes(filter_=['HOXC11', 'HOXC13', 'HNRNPA1', 'CBX5'], label_aligns={'HOXC11':'left', 'HOXC13':'right', 'CBX5':'right'}, fontsize=40)
vis.plot_signal('RNA-Seq', 'data_HOXC13/235990-RN01.Aligned.sortedByCoord.out.bam.bw', label_size=45, data_range_size=45, max_value=3500, color='#E31A1C')
vis.plot_signal('ATAC-Seq', 'data_HOXC13/1459011-34_ATAC.pval.signal.bigwig', label_size=45, data_range_size=45, max_value=100, color='#6A3D9A')
vis.plot_arcs(lw=4, cutoff='top', gene_filter=['HOXC11', 'HOXC13'], arc_color='#666666') # IKZF3-related neo-loops
vis.plot_chromosome_bar(name_size=50, coord_size=45, color_by_order=['#1F78B4','#33A02C'])

vis.fig.set_facecolor('white')
