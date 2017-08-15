import numpy as np

def resize_matrix(a,coef=10,oper=np.average):
	l=len(a)
	res = []
	for i in range(0,l/coef):
		res.append([])
		for j in range(0,l/coef):
			i_border=(i+1)*(coef)
			j_border=(j+1)*(coef)
			if i_border > l: i_border=l
			if j_border > l: j_border=l
#			print i*coef,i_border,j*coef,j_border
#			print a[i*coef:i_border,j*coef:j_border]
#			print np.average(a[i*coef:i_border,j*coef:j_border])
			elem = oper(a[i*coef:i_border,j*coef:j_border])
			res[i].append(elem)
	return np.array(res)
	
def get_idx_of_chrms(g): #returns array of arrays, containig indexes of chrms in genome matrix
	res = []
	for i in xrange(g.chrmCount):
		res.append (range(g.chrmStartsBinCont[i],g.chrmEndsBinCont[i]))

	return np.array(res)
