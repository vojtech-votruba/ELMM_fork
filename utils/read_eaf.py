#!/usr/bin/python3
import numpy as np

scalar_flag = 1001

vector_flag = 1002
   

class eaf_reader:
  
    def __init__(self, base_name, bounding_box=None):
    
        self.base_name = base_name
        
        self.bounding_box = bounding_box
        if self.bounding_box is not None:
            self.bounding_box = np.asarray(self.bounding_box)
        

        self.efh_name = base_name + ".efh"

        efh = open(self.efh_name, "rb")

        one = np.fromfile(efh, dtype="int32", count=1)[0]

        if (one != 1):
            print("Error, the data has the the opposite endianness!")
            raise ValueError
            
        
        float_size = np.fromfile(efh, dtype="int32", count=1)[0]
        
        if (float_size == 4):
            self.dtype = "float32"
        elif (float_size == 8):
            self.dtype = "float64"
        else:
            print("Unknown float size in " + efh_name)
            raise ValueError
          
        self.nxyz = np.fromfile(efh, dtype="int32", count=3).tolist()
        self.orig_nxyz = self.nxyz

        
        self.xs = np.fromfile(efh, dtype = self.dtype, count = self.nxyz[0])
        self.ys = np.fromfile(efh, dtype = self.dtype, count = self.nxyz[1])
        self.zs = np.fromfile(efh, dtype = self.dtype, count = self.nxyz[2])
        
        if self.bounding_box is not None:
            for i in range(self.xs.size):
                self.mini = self.xs.size
                if self.xs[i] >= self.bounding_box[0]:
                   self.mini = i
                   break
            for j in range(self.ys.size):
                self.minj = self.ys.size
                if self.ys[j] >= self.bounding_box[2]:
                   self.minj = j
                   break
            for k in range(self.zs.size):
                self.mink = self.zs.size
                if self.zs[k] >= self.bounding_box[4]:
                   self.mink = k
                   break
                    
            for i in range(self.xs.size-1,-1,-1):
                self.maxi = 0
                if self.xs[i] <= self.bounding_box[1]:
                   self.maxi = i+1
                   break
            for j in range(self.ys.size-1,-1,-1):
                self.maxj = 0
                if self.ys[j] <= self.bounding_box[3]:
                   self.maxj = j+1
                   break
            for k in range(self.zs.size-1,-1,-1):
                self.maxk = 0
                if self.zs[k] <= self.bounding_box[5]:
                   self.maxk = k+1
                   break
                 
            self.nxyz = [self.maxi-self.mini,
                         self.maxj-self.minj,
                         self.maxk-self.mink]
            self.xs = self.xs[self.mini:self.maxi]
            self.ys = self.ys[self.minj:self.maxj]
            self.zs = self.zs[self.mink:self.maxk]
                    
        else:
            self.mini = 0
            self.maxi = self.xs.size #intentionally +1
            self.minj = 0
            self.maxj = self.ys.size
            self.mink = 0
            self.maxk = self.zs.size
        
        self.arrays = []
        #Now for simplicity just assume that we know that the file contains velocity vectors only.
        #Therefore we skip reading the data description.
        self.arrays = self.arrays + [{"name":"u", "vector":True}]
        
        self.transformed_xy = False
        
        
        
    def read_frame(self, n_frame):
        eaf = open(self.base_name + "-" + str(n_frame) + ".eaf", "rb")
        
        data = []
        
        for array in self.arrays:
            if array["vector"]:
              array_data = np.fromfile(eaf, dtype=self.dtype).reshape([3]+self.orig_nxyz, order='F')
              
              array_data = array_data[:,
                                      self.mini:self.maxi,
                                      self.minj:self.maxj,
                                      self.mink:self.maxk]
              
              if self.transformed_xy:
                array_data = array_data.swapaxes(1,2)
                array_data[[0,1],:,:,:] = array_data[[1,0],:,:,:]            

            else:
              array_data = np.fromfile(eaf, dtype=self.dtype).reshape(self.orig_nxyz, order='F')
              
              array_data = array_data[self.mini:self.maxi,
                                      self.minj:self.maxj,
                                      self.mink:self.maxk]
              
              if self.transformed_xy:
                array_data = array_data.swapaxes(0,1)
              
            data = data + [{"name" : array["name"], "data" : array_data}]
            
        return data
   
    #returns the x-xoordinate of the ith index
    def x(self,i):
      return self.xs[i]
      
    #returns the y-xoordinate of the ith index
    def y(self,i):
      return self.ys[i]
      
    #returns the z-xoordinate of the ith index
    def z(self,i):
      return self.zs[i]
    
    #converts the vector components u,v,w in data to one long 1D sequence
    #for certain specialized purposes, ignore if not needed
    def step_data_matrix(self, data, name, swap=False):
        n = np.product(self.nxyz)
        for i in range(3*n):
          if data[0]["name"] == name:
            u = data[0]["data"]
            if swap is not None:
                u = u.swapaxes(swap[0],swap[1])
            res = np.ndarray([3*n])
            res[0:n] = u[0,:,:,:].reshape([n])
            res[n:2*n] = u[1,:,:,:].reshape([n])
            res[2*n:3*n] = u[2,:,:,:].reshape([n])
            return res
        return None    
     
    #returns a big array with a row of step_data_matrix times each time step
    def steps_data_matrix(self, name, first, last, step=1, swap=None):
        n = np.product(self.nxyz)
        nsteps = (last - first + 1) // step
        res = np.ndarray([3*n, nsteps])
        
        k = 0
        for tstep in range(first, last+1, step):
          print(tstep)
          data = self.read_frame(tstep)
          res[:,k] = self.step_data_matrix(data, name, swap)
          k += 1
        return res
          
          
    #switch x and y (and u and v)
    #useful when simulationg a channel in the y direction
    def transform_xy(self):
        self.xs, self.ys = self.ys, self.xs
        self.nxyz = np.asarray([self.nxyz[1],self.nxyz[0],self.nxyz[2]])
        self.transformed_xy = True
        
        #if self.bounding_box is not None:
            #self.bounding_box[[1,0]] = self.bounding_box[[0,1]]
            
            #for i in range(self.xs.size):
                #self.mini = self.xs.size
                #if self.xs[i] >= self.bounding_box[0]:
                   #self.mini = i
                   #break
            #for j in range(self.ys.size):
                #self.minj = self.ys.size
                #if self.ys[j] >= self.bounding_box[2]:
                   #self.minj = j
                   #break
            #for k in range(self.zs.size):
                #self.mink = self.zs.size
                #if self.zs[k] >= self.bounding_box[4]:
                   #self.mink = k
                   #break
                    
            #for i in range(self.xs.size-1,-1,-1):
                #self.maxi = 0
                #if self.xs[i] <= self.bounding_box[1]:
                   #self.maxi = i+1
                   #break
            #for j in range(self.ys.size-1,-1,-1):
                #self.maxj = 0
                #if self.ys[j] <= self.bounding_box[3]:
                   #self.maxj = j+1
                   #break
            #for k in range(self.zs.size-1,-1,-1):
                #self.maxk = 0
                #if self.zs[k] <= self.bounding_box[5]:
                   #self.maxk = k+1
                   #break
                    
        #else:
            #self.mini = 0
            #self.maxi = self.xs.size #intentionally +1
            #self.minj = 0
            #self.maxj = self.ys.size
            #self.mink = 0
            #self.maxk = self.zs.size
       
      
                           
if __name__ == '__main__':
   name = "frame-test"
   
   reader = eaf_reader(name)
   
   print(reader.x(0))
   print(reader.y(0))
   print(reader.z(0))
   
   data = reader.read_frame(0)
   
   print(data[0]["name"])
   
   u = data[0]["data"]
   
   print(np.shape(u))
   
   #first index for vectors: 0..x component, 1..y component, 2..z component
   
   #then x index, y index, z index
   print(u[0,0,0,0])
   print(u[1,1,0,0])
   print(u[2,0,0,1])
   
   
   
   
