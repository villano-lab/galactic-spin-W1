import re

"""various functions to help deal with text data in Python """
def getXYdata_wXYerr(infile=None):

  #format expected is (whitespace separated):
  #x	y	yerr

  #can handle comment lines where the line's first non-whitespace character is #

  #open the file return a dictionary with se,sw,nw,ne elements
  f = open(infile)

  #make a list for vector identifier 
  #first two are x-y of histogram-type step function
  #second two are a sort-of smooth curve to represent the function
  vecs = ['xx','yy','ex','ey']

  #make a dictionary to store the pulses
  funcs = {}

  #read file N times, is this efficient?
  regex=re.compile(r'^\s*#.+')
  #[print(regex.search(x)) for x in f.readlines()]
  funcs[vecs[0]] = [x.split()[0] for x in f.readlines() if regex.search(x) is None]
  f.seek(0)
  funcs[vecs[1]] = [x.split()[1] for x in f.readlines() if regex.search(x) is None]
  f.seek(0)
  funcs[vecs[2]] = [x.split()[2] for x in f.readlines() if regex.search(x) is None]
  f.seek(0)
  funcs[vecs[3]] = [x.split()[3] for x in f.readlines() if regex.search(x) is None]

  f.close()

  #convert to floats
  funcs[vecs[0]] = [float(i) for i in funcs[vecs[0]]]
  funcs[vecs[1]] = [float(i) for i in funcs[vecs[1]]]
  funcs[vecs[2]] = [float(i) for i in funcs[vecs[2]]]
  funcs[vecs[3]] = [float(i) for i in funcs[vecs[3]]]

  return funcs 

def getXYdata_wYerr(infile=None):

  #format expected is (whitespace separated):
  #x	y	yerr

  #can handle comment lines where the line's first non-whitespace character is #

  #open the file return a dictionary with se,sw,nw,ne elements
  f = open(infile)

  #make a list for vector identifier 
  #first two are x-y of histogram-type step function
  #second two are a sort-of smooth curve to represent the function
  vecs = ['xx','yy','ey']

  #make a dictionary to store the pulses
  funcs = {}

  #read file N times, is this efficient?
  regex=re.compile(r'^\s*#.+')
  #[print(regex.search(x)) for x in f.readlines()]
  funcs[vecs[0]] = [x.split()[0] for x in f.readlines() if regex.search(x) is None]
  f.seek(0)
  funcs[vecs[1]] = [x.split()[1] for x in f.readlines() if regex.search(x) is None]
  f.seek(0)
  funcs[vecs[2]] = [x.split()[2] for x in f.readlines() if regex.search(x) is None]

  f.close()

  #convert to floats
  funcs[vecs[0]] = [float(i) for i in funcs[vecs[0]]]
  funcs[vecs[1]] = [float(i) for i in funcs[vecs[1]]]
  funcs[vecs[2]] = [float(i) for i in funcs[vecs[2]]]

  return funcs 
def getXYdata(infile=None):

  #format expected is (whitespace separated):
  #x	y

  #can handle comment lines where the line's first non-whitespace character is #

  #open the file return a dictionary with se,sw,nw,ne elements
  f = open(infile)

  #make a list for vector identifier 
  #first two are x-y of histogram-type step function
  #second two are a sort-of smooth curve to represent the function
  vecs = ['xx','yy']

  #make a dictionary to store the pulses
  funcs = {}

  #read file N times, is this efficient?
  regex=re.compile(r'^\s*#.+')
  #[print(regex.search(x)) for x in f.readlines()]
  funcs[vecs[0]] = [x.split()[0] for x in f.readlines() if regex.search(x) is None]
  f.seek(0)
  funcs[vecs[1]] = [x.split()[1] for x in f.readlines() if regex.search(x) is None]

  f.close()

  #convert to floats
  funcs[vecs[0]] = [float(i) for i in funcs[vecs[0]]]
  funcs[vecs[1]] = [float(i) for i in funcs[vecs[1]]]

  return funcs 
def getXYZdata(infile=None):

  #format expected is (whitespace separated):
  #x	y   z

  #can handle comment lines where the line's first non-whitespace character is #

  #open the file return a dictionary with se,sw,nw,ne elements
  f = open(infile)

  #make a list for vector identifier 
  #first two are x-y of histogram-type step function
  #second two are a sort-of smooth curve to represent the function
  vecs = ['xx','yy','zz']

  #make a dictionary to store the pulses
  funcs = {}

  #read file N times, is this efficient?
  regex=re.compile(r'^\s*#.+')
  #[print(regex.search(x)) for x in f.readlines()]
  funcs[vecs[0]] = [x.split()[0] for x in f.readlines() if regex.search(x) is None]
  f.seek(0)
  funcs[vecs[1]] = [x.split()[1] for x in f.readlines() if regex.search(x) is None]
  f.seek(0)
  funcs[vecs[2]] = [x.split()[2] for x in f.readlines() if regex.search(x) is None]

  f.close()

  #convert to floats
  funcs[vecs[0]] = [float(i) for i in funcs[vecs[0]]]
  funcs[vecs[1]] = [float(i) for i in funcs[vecs[1]]]
  funcs[vecs[2]] = [float(i) for i in funcs[vecs[2]]]

  return funcs 
def getXYZABCdata(infile=None):

  #format expected is (whitespace separated):
  #x	y   z	a   b   c

  #can handle comment lines where the line's first non-whitespace character is #

  #open the file return a dictionary with se,sw,nw,ne elements
  f = open(infile)

  #make a list for vector identifier 
  #first two are x-y of histogram-type step function
  #second two are a sort-of smooth curve to represent the function
  vecs = ['xx','yy','zz','aa','bb','cc']

  #make a dictionary to store the pulses
  funcs = {}

  #read file N times, is this efficient?
  regex=re.compile(r'^\s*#.+')
  #[print(regex.search(x)) for x in f.readlines()]
  funcs[vecs[0]] = [x.split()[0] for x in f.readlines() if regex.search(x) is None]
  f.seek(0)
  funcs[vecs[1]] = [x.split()[1] for x in f.readlines() if regex.search(x) is None]
  f.seek(0)
  funcs[vecs[2]] = [x.split()[2] for x in f.readlines() if regex.search(x) is None]
  f.seek(0)
  funcs[vecs[3]] = [x.split()[3] for x in f.readlines() if regex.search(x) is None]
  f.seek(0)
  funcs[vecs[4]] = [x.split()[4] for x in f.readlines() if regex.search(x) is None]
  f.seek(0)
  funcs[vecs[5]] = [x.split()[5] for x in f.readlines() if regex.search(x) is None]

  f.close()

  #convert to floats
  funcs[vecs[0]] = [float(i) for i in funcs[vecs[0]]]
  funcs[vecs[1]] = [float(i) for i in funcs[vecs[1]]]
  funcs[vecs[2]] = [float(i) for i in funcs[vecs[2]]]
  funcs[vecs[3]] = [float(i) for i in funcs[vecs[3]]]
  funcs[vecs[4]] = [float(i) for i in funcs[vecs[4]]]
  funcs[vecs[5]] = [float(i) for i in funcs[vecs[5]]]

  return funcs 
