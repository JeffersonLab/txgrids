try:
    import _pickle as cPickle
except:
    import cPickle
import sys, os
import zlib

#--physics 

def convert_lum(lum):
    one=0.3893793721  #--GeV2 mbarn from PDG
    lum,units=lum.split(':')
    lum=float(lum)
    units=units.strip()
    if units=='fb-1':   return lum*one*1e12
    else: sys.exit('units not convertible!')

#--aux 

def checkdir(path):
  if not os.path.exists(path): 
    os.makedirs(path)

def save(data,name):  
  #compressed=zlib.compress(cPickle.dumps(data))
  #f=open(name,"wb")
  #f.writelines(compressed)
  #f.close()
  cPickle.dumps(data)
  compressed=zlib.compress(cPickle.dumps(data))
  f=open(name,"wb")
  f.write(compressed)
  f.close()

def load(name): 
  compressed=open(name,"rb").read()
  data=cPickle.loads(zlib.decompress(compressed))
  return data

def lprint(msg):
  sys.stdout.write('\r')
  sys.stdout.write(msg)
  sys.stdout.flush()



