from sys import argv
from collections import Counter, defaultdict
import numpy as np
import matplotlib.pyplot as plt
import re
import fileinput
import math
import time

################################################################
def voc():
	dict = defaultdict()
	dict.default_factory = lambda:len(dict)
	return dict

def get_ids(C, dict):
	yield[dict[w] for w in C]

def cosSimilarity(a,b):
	#if (len(a) != len(b)): return 666
	m1 = 0
	m2 = 0
	for i in xrange(a.shape[0]):
		m1 = m1 + a[i]*a[i]
		m2 = m2 + b[i]*b[i]
	m1 = math.sqrt(m1)
	m2 = math.sqrt(m2)
	return np.dot(a,b)/(m1*m2)

def tanhDist(a,b):
	return np.tanh(np.dot(a,b))

def probContexto(W, V, obj = 0, ctx=0):    #ctx=context
	n = V.shape[1]
	acum = 0
	for i in xrange(n):
		if i == obj: continue
		aux = tanhDist(W[i,:], ctx)
		acum += math.exp(-aux)
	aux = math.exp(-1*tanhDist(W[obj,:], ctx))
	return aux/acum

def negativeSampling(W,V, obj=0, ctx=0, member=1): 
	hac = ctx * probContexto(W, V, obj, ctx)
	if member:
		return 1-hac
	else:
		return -hac
	
def replaceV(V, v, etha=0.2, obj=0):
	V[:,obj] = V[:,obj] - etha*np.matrix((v))
	return V

def actualizaW(W, V, etha = 0.2, contexto = 0):
	#print W, contexto
	#print '*********************'
	hac = np.zeros(2)
	for i in xrange(V.shape[1]):
		for j in xrange(V.shape[0]):
			hac[j] = hac[j]+V[j,i]
	hac = hac/V.shape[1]
	hac = hac*etha
	for i in xrange(len(contexto)):
		pos = W_v[contexto[i]]
		#print contexto[i], pos, W[pos, :]
		W[pos, :] = W[pos,:]- etha*np.matrix((hac))
		#print W
	return W
################################################################

np.set_printoptions(threshold='nan')
word = re.compile(r'[\w]+', re.UNICODE) 	  #expresion regular para encontrar palabras
count = 0                           		 #numero de palabras en el documento
W_f = Counter()               			     #cuantas palabras de cada tipo, diccionario

for line in fileinput.input():
	aux = re.findall(word, line.decode('utf8'))
	for i in xrange(len(aux)):
		W_f[aux[i]] += 1
		count += 1

W_v = voc()
list(get_ids(W_f.keys(), W_v))[0]                #diccionario de palbras

tipos = len(W_v)                                 #cantidad de palabras diferentes (tipos)
dimension = 2                                    #referente al espacio vectorial
#nn = np.array(range(tipos))                      #etiqueas arbitrarias para graficar 
#nn = 'bebe el caballo pezcado saluda gato perro corre brinca come raton'
nn = 'el gato perro corre come raton'
nn = nn.split()
print W_v

W = 10*np.random.rand(tipos,dimension)          #Matrices iniciales, aleatorias
V = 10*np.random.rand(dimension,tipos)

plt.ion()                                       #grafica la distribucion de la matriz W inicial (aleatoria)
fig, ax = plt.subplots()
ax.scatter(W[:,0], W[:,1])
for i, txt in enumerate(nn):
	ax.annotate(txt, (W[i,0], W[i,1]))
plt.pause(0.0001)
time.sleep(1)

mem = 1           #booleano, para determinar si forma parte del contexto la palabras
iteracion = 20    #numero de iteraciones que correran

for k in xrange(iteracion):
	for words in W_v.keys():
		for line in fileinput.input():
			hcomp = 0
			aux = re.findall(word, line.decode('utf8'))
			if len(aux) == 0: continue
			for i in xrange(len(aux)):
				hcomp += V[:, W_v[aux[i]]]
			hcomp = hcomp/len(aux)
			if word in aux:
				mem = 1
			else:
				mem = 0
			for i in xrange(len(aux)):
				aux2 = negativeSampling(W, V, obj=W_v[words], ctx=hcomp, member=mem)
				replaceV(V, aux2, obj=W_v[aux[i]])
			actualizaW(W,V,contexto = aux)

plt.ion()
fig, ax = plt.subplots()
ax.scatter(W[:,0], W[:,1])
for i, txt in enumerate(nn):
	ax.annotate(txt, (W[i,0], W[i,1]))
plt.pause(0.0001)
time.sleep(5)

print W_v