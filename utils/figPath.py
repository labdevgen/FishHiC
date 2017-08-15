import sys,os

figure_path = sys.argv[0]+"_figs/"
if not os.path.exists(figure_path):
	os.makedirs(figure_path)