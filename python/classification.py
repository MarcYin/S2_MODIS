import sys
sys.path.insert(0,'./python')
from readSent import *
import numpy 
from numpy import clip, where
from scipy.ndimage.morphology import *
import xml.etree.cElementTree as ET
#import pathos.multiprocessing as mp


class classification(object):
    
    def __init__(self, img = None, fhead = 'data/50SMG20164100', bands = (2,3,4,8,11,12,13), bounds = None):
        self.img = img
        self.bounds = bounds 
        self.fhead = fhead
        self.bands = bands
        self.t1_max = 0.25
        self.t1_min = 0.07
        self.t2_max = 0.2
        self.t2_min = -0.1
        self.t3_max = 0.42
        self.t3_min = 0.2
        self.t4_max = 0.35
        self.t4_min = 0.15
        self.t5_max = 0.22
        self.t5_min = 0.18
        self.t6_max = 0.95
        self.t6_min = 0.85
        self.t7_eb12 = 0.35
        self.t7_sp = 0.73 #(0.9^5 5 means 5 processes)
        self.t8_max = 0.2
        self.t8_min = 0.02
        self.t9_max = 1.6
        self.t9_min = 1.
        self.t10_max = 0.9
        self.t10_min = 0.6
        self.t11_1 = 0.034
        self.t11_2 = 0.2
        self.t12_max = 0.21
        self.t12_min = 0.07
        self.t13_max = 1.1
        self.t13_min = 0.9
        self.t14_max = 6.
        self.t14_min = 3.
        self.c_pt = 0.015
        '''self.b2 = img['B02']
        self.b4 = img['B04']
        self.b3 = img['B03']
        self.b8 = img['B08']
        self.b11 = img['B11']
        self.b12 = img['B12']
        self.b8a = img['B8A']
        '''   
    def Ratio(self, data, minm, maxm): # get the posibility of the values within specific bound
        _ratio = (data-minm)/(maxm-minm)
        return _ratio
    
    def Posb(self, p1, p2): # get the combined possibility of two array
        _posb = p1*p2
        return _posb
    
    def CombinM(self, m1, m2,t = True): # combination of two mask either & or | operator
        if t:
            m = m1 & m2
        else:
            m = m1 | m2  
        return m

    def Getimg(self,bands, fhead, bounds): # read in the img data
        #print 'get some imgs'
        print fhead, bands
        img = readfile(bands, fhead, bounds)
        self.b2 = img['B02']
        self.b4 = img['B04']
        self.b3 = img['B03']
        self.b8 = img['B08']
        self.b11 = img['B11']
        self.b12 = img['B12']
        self.b8a = img['B8A'] 
        del img
    
    def ScaleExtent(self, data): # used for unifine different array, 

        shape = self.b2.shape
        re = int(shape[0]/(data.shape[0]))
        return numpy.repeat(numpy.repeat(data, re, axis = 1), re, axis =0)

    def M_P_func(self, data, minm, maxm, nm = 1):
        ca = clip(data, minm, maxm)
        r = self.Ratio(ca, minm, maxm)
        if nm == 1:
            m = ca > minm
        elif nm == 2:
            m = ca < minm
        elif nm == 3:
            m = ca > maxm
        elif nm == 4:
            m = ca < maxm
        else:
            pass
        return (r, m)
    
    def _del(self, v):
        for i in v:
            del i

    def _p1(self): # process 1 using the band 4 with thresh 1 (t1_max and t1_min)
        
        self.p1p, self.p1m = self.M_P_func(self.b4, self.t1_min, self.t1_max)
        return self.p1p,self.p1m   


    def _p1(self): # process 1 using the band 4 with thresh 1 (t1_max and t1_min)
        
        self.p1p, self.p1m = self.M_P_func(self.b4, self.t1_min, self.t1_max)
        return self.p1p,self.p1m   

    def _p2(self): # precess 2 using the NDSI with thresh 2
        
        self.eb11 = self.ScaleExtent(self.b11)
        self.NDSI = (self.b3 - self.eb11)/(self.eb11 + self.b3)
        self.p2p, self.p2m = self.M_P_func(self.NDSI, self.t2_min, self.t2_max)
        return self.p2p,self.p2m
        self._del(self.b11)

    # snow detection and get the cloud mask  
    
    def _p3(self): # NDSI snow detection
        
        self.p3p, self.p3m = self.M_P_func(self.NDSI, self.t3_min, self.t3_max)
        return self.p3p, self.p3m

    
    def _p4(self): # b8 snow
        
        self.p4p, self.p4m = self.M_P_func(self.b8, self.t4_min, self.t4_max)
        return self.p4p, self.p4m


    def _p5(self): # b2 snow
        
        self.p5p, self.p5m = self.M_P_func(self.b2, self.t5_min, self.t5_max)
        return self.p5p, self.p5m
    
    
    def _p6(self): # b2/b4 snow
        
        b2_b4 = self.b2/self.b4
        
        self.p6p, self.p6m = self.M_P_func(b2_b4, self.t6_min, self.t6_max)
        return self.p6p, self.p6m

    
    
    
    def _p7(self): # dilation for snow adjacent pixels 
        
        self.sm = numpy.ones(self.b2.shape).astype(bool)
       
        for i in [self.p3m, self.p4m, self.p5m, self.p6m]:
            self.sm = self.CombinM(i, self.sm)
        
        self.eb12 = self.ScaleExtent(self.b12)
        self._del(self.b12)
        b12m = self.eb12 < self.t7_eb12
        
        struct = iterate_structure(generate_binary_structure(2,1), 3)
        snow_mask_dil = binary_dilation(self.sm, struct)
        bz = snow_mask_dil -  self.sm
        bs = bz & b12m
        #print bs
        self.sm = self.CombinM(bs, self.sm, t=False)
        
        self.sp = numpy.ones(self.b2.shape)
        
        for i in [self.p3p, self.p4p, self.p5p, self.p6p]:
            self.sp = self.Posb(i,self.sp)
            
        
        S_p_m = self.sp > self.t7_sp #use the finnal possibility to thershold (0.73 for his) the snow mask
        self.sm = self.CombinM(self.sm, S_p_m)

    def _p8(self):
        
        self.NDVI = (self.b8 - self.b4)/(self.b8 + self.b4)
        p8p, self.p8m = self.M_P_func(self.NDVI, self.t8_min, self.t8_max, nm = 4) 
        self.p8p =  p8p*-1 + 1
        self.vm = ~ self.p8m

    def _p9(self):
        
        b8_b3 = self.b8/self.b3
        p9p, self.p9m = self.M_P_func(b8_b3, self.t9_min, self.t9_max, nm=4)
        self.p9p = 1. - p9p
        return self.p9p, self.p9m

    def _p10(self):

        b2_b11 = self.b2/self.eb11
        
        b2_FT = clip(-0.4*b2_b11 + 0.46, 0.15, 0.32)        
        
        self.bsm = (self.b2 < b2_FT) & (b2_b11 < 0.6) #bare soil
        self.p10p, p10m = self.M_P_func(b2_b11, self.t10_min, self.t10_max)
        self.p10m = p10m & ~(self.b2 < b2_FT)

    # water bodies
    def _p11(self):
        
        self.eb8a = self.ScaleExtent(self.b8a)
        d_b2_b4 = self.b2-self.b4
        self.wm = (d_b2_b4 > self.t11_1) & (self.eb8a < self.b4) & (self.b2 < self.t11_2)
        self.p11m = ~self.wm 
        self.p11p = 1.
        


    
    def _p12(self):
        
        b2_b11 = self.b2/self.eb11
        B12_FT = clip(b2_b11*0.1 - 0.09, self.t12_min, self.t12_max)
        self.wm2 = where((b2_b11 > 4.) & (self.eb12 < B12_FT) & (self.eb8a < self.b4) & (self.b2 < 0.2), True, False)
        R15_AMB = (b2_b11 < 4.) & (b2_b11 >= 2.) & (self.eb12 < B12_FT) & (self.eb8a < self.b4) & (self.b2 < 0.2)
        self.wm3 = R15_AMB | self.wm | self.wm2
        self.p12m = ~self.wm3
        self.p12p = 1.

    # sand and stone in desert
    def _p13(self):
        
        self.eb8a = self.ScaleExtent(self.b8a)
        b8a_b11 = self.eb8a/self.eb11
        b2_FT = clip(b8a_b11*(-0.25)+0.475, 0.16, 0.35)
        self.Dssm = (self.b2 < b2_FT) & (b8a_b11 < 0.9) & (self.b2 < 0.8)
        self.p13p, p13m = self.M_P_func(b8a_b11, self.t13_min, self.t13_max)
        self.p13m = p13m & ~(self.b2 < b2_FT)
        return self.p13p, self.p13m, self.Dssm 

    def _p14(self):
        
        b4_b11 = self.b4/self.eb11
        p14p, self.p14m = self.M_P_func(b4_b11, self.t14_min, self.t14_max, nm = 4)
        self.p14p = 1. - p14p
        return self.p14p,self.p14m

    def Get_cm_p(self):
        print '\nHere is the cloud probability calculation!\n'
        if self.img == None:
            
            self.Getimg(self.bands, self.fhead, self.bounds)
        else:
            pass
        self._p1()
        print "{0:.0f}%".format(1./14 * 100)
        self._p2()
        print "{0:.0f}%".format(2./14 * 100)
        self._p3()
        print "{0:.0f}%".format(3./14 * 100)
        self._p4()
        print "{0:.0f}%".format(4./14 * 100)
        self._p5()
        print "{0:.0f}%".format(5./14 * 100)
        self._p6()
        print "{0:.0f}%".format(6./14 * 100)
        self._p7()
        print "{0:.0f}%".format(7./14 * 100)
        self._p8()
        print "{0:.0f}%".format(8./14 * 100)
        self._p9()
        print "{0:.0f}%".format(9./14 * 100)
        self._p10()
        print "{0:.0f}%".format(10./14 * 100)
        self._p11()
        print "{0:.0f}%".format(11./14 * 100)
        self._p12()
        print "{0:.0f}%".format(12./14 * 100)
        self._p13()
        print "{0:.0f}%".format(13./14 * 100)
        self._p14()
        print "{0:.0f}%".format(14./14 * 100-1.111)
       
        self.cm = True
        for i in [self.p1m,self.p2m,~self.sm,self.p8m, self.p9m, self.p10m, self.p11m, self.p12m, self.p13m, self.p14m]:
            self.cm = self.CombinM(i, self.cm)
   
        
        self.cp = 1.
        for i in [self.p1p, self.p2p, self.p8p, self.p9p, self.p10p, self.p11p, self.p12p, self.p13p, self.p14p]:
            self.cp = self.Posb(i,self.cp)
        
        c_ptm = self.cp > self.c_pt
        self.fcm = self.cm & c_ptm
        
        v = [self.p8m, self.p9m, self.p10m, self.p11m, self.p12m, self.p13m, self.p14m, 
                  self.p8p, self.p9p, self.p10p, self.p11p, self.p12p, self.p13p, self.p14p]
        self._del(v)
        
        print 'Done!!!'


    
    
'''
{k.strip(): float(v) for k, v in a}
{'self.c_pt': 0.015,
 'self.t10_max': 0.9,
 'self.t10_min': 0.6,
 'self.t11_1': 0.034,
 'self.t11_2': 0.2,
 'self.t12_max': 0.21,
 'self.t12_min': 0.07,
 'self.t13_max': 1.1,
 'self.t13_min': 0.9,
 'self.t14_max': 6.0,
 'self.t14_min': 3.0,
 'self.t1_max': 0.25,
 'self.t1_min': 0.07,
 'self.t2_max': 0.2,
 'self.t2_min': -0.1,
 'self.t3_max': 0.42,
 'self.t3_min': 0.2,
 'self.t4_max': 0.35,
 'self.t4_min': 0.15,
 'self.t5_max': 0.22,
 'self.t5_min': 0.18,
 'self.t6_max': 0.95,
 'self.t6_min': 0.85,
 'self.t7_eb12': 0.35,
 'self.t7_sp': 0.73,
 'self.t8_max': 0.2,
 'self.t8_min': 0.02,
 'self.t9_max': 1.6,
 'self.t9_min': 1.0}
 
       


    


'''
"""
class classification(object):
    
    def __init__(self, fhead = 'data/50SMG20164100', bands = (2,3,4,8,11,12,13), bounds = None):
        self.bounds = bounds 
        self.fhead = fhead
        self.bands = bands
        self.t1_max = 0.25
        self.t1_min = 0.07
        self.t2_max = 0.2
        self.t2_min = -0.1
        self.t3_max = 0.42
        self.t3_min = 0.2
        self.t4_max = 0.35
        self.t4_min = 0.15
        self.t5_max = 0.22
        self.t5_min = 0.18
        self.t6_max = 0.95
        self.t6_min = 0.85
        self.t7_eb12 = 0.35
        self.t7_sp = 0.73 #(0.9^5 5 means 5 processes)
        self.t8_max = 0.2
        self.t8_min = 0.02
        self.t9_max = 1.6
        self.t9_min = 1.
        #self.t10_max = 0.32
        #self.t10_min = 0.15
        self.t10_max = 0.9
        self.t10_min = 0.6
        self.t11_1 = 0.034
        self.t11_2 = 0.2
        self.t12_max = 0.21
        self.t12_min = 0.07
        self.t13_max = 1.1
        self.t13_min = 0.9
        self.t14_max = 6.
        self.t14_min = 3.
        self.c_pt = 0.015
        #self.t15_max = 
        #self.t15_min = 0.07
        #self.t16_max = 0.35
        #self.t16_min = 0.16
        self.cm = []
        self.sm = []
        self.wm = []
        self.bsm= []
        self.vm = []
        self.tempM = []
        self.Dssm = []
        '''
        self.b2 = img['B02']
        self.b4 = img['B04']
        self.b3 = img['B03']
        self.b8 = img['B08']
        self.b11 = img['B11']
        self.b12 = img['B12']
        self.b8a = img['B8A']
        ''' 
        
        
        
        
    def Ratio(self, data, minm, maxm): # get the posibility of the values within specific bound
        _ratio = (data-minm)/(maxm-minm)
        return _ratio
    
    def Posb(self, p1, p2): # get the combined possibility of two array
        _posb = p1*p2
        return _posb
    
    def CombinM(self, m1, m2,t = True): # combination of two mask either & or | operator
        if t:
            m = m1 & m2
        else:
            m = m1 | m2
        
        return m
    
    def Getimg(self,bands, fhead, bounds): # read in the img data
        print 'get some imgs'
        print fhead, bands
        img = readfile(bands, fhead, bounds)
 
      
        self.b2 = img['B02']
        self.b4 = img['B04']
        self.b3 = img['B03']
        self.b8 = img['B08']
        self.b11 = img['B11']
        self.b12 = img['B12']
        self.b8a = img['B8A']
        
        del img
    
    def ScaleExtent(self, data): # used for unifine different array, 


        shape = self.b2.shape
        re = int(shape[0]/(data.shape[0]))
        return numpy.repeat(numpy.repeat(data, re, axis = 1), re, axis =0)


    def M_P_func(self, data, minm, maxm, nm = 1):
        ca = clip(data, minm, maxm)
        r = self.Ratio(ca, minm, maxm)
        if nm == 1:
            m = ca > minm
        elif nm == 2:
            m = ca < minm
        elif nm == 3:
            m = ca > maxm
        elif nm == 4:
            m = ca < maxm
        else:
            pass
        return (r, m)
    
    def _del(self, v):
        for i in v:
            del i
    
    def _p1(self): # process 1 using the band 4 with thresh 1 (t1_max and t1_min)
        
        self.p1p, self.p1m = self.M_P_func(self.b4, self.t1_min, self.t1_max)
        #ca = clip(self.b4, self.t1_min, self.t1_max)
        #self.Ratio(ca, self.t1_min, self.t1_max)
        #self.tempP = self.ratio
        #self.tempcm = ca>self.t1_min
    
    def _p2(self): # precess 2 using the NDSI with thresh 2
        
        self.eb11 = self.ScaleExtent(self.b11)
        self.NDSI = (self.b3 - self.eb11)/(self.eb11 + self.b3)
        self.p2p, self.p2m = self.M_P_func(self.NDSI, self.t2_min, self.t2_max)
        #ca = clip(self.NDSI, sefl.t2_min, self.t2_max)
        #self.Ratio(ca, self.t2_min, self.t2_max)
        self._del(self.b11)
    
        # combination of previous 2 stadges
        #self.p1_p2_P = self.Posb(self.p1p, self.p2p)
        #self.p1_p2_M = self.CombinM(self.p2m, self.p1m)
        #self.cm = sefl.combinm
        #self.tempP = self.ratio
        #self._del([self.p1p, self.p2p, self.p1m, self.p2m]) # delet variables
    
    
    
    # snow detection and get the cloud mask  
    
    def _p3(self): # NDSI snow detection
        
        self.p3p, self.p3m = self.M_P_func(self.NDSI, self.t3_min, self.t3_max)
        #ca = clip(self.NDSI, sefl.t3_min, self.t3_max)
        #self.Ratio(ca, self.t3_min, self.t3_max)
        #self.sratio = self.ratio
        #self.tempsm = ca > self.t3_min
    
    def _p4(self): # b8 snow
        
        self.p4p, self.p4m = self.M_P_func(self.b8, self.t4_min, self.t4_max)
    
    def _p5(self): # b2 snow
        
        self.p5p, self.p5m = self.M_P_func(self.b2, self.t5_min, self.t5_max)
    
    def _p6(self): # b2/b4 snow
        
        b2_b4 = self.b2/self.b4
        
        self.p6p, self.p6m = self.M_P_func(b2_b4, self.t6_min, self.t6_max)
    
    
    
    
    def _p7(self): # dilation for snow adjacent pixels 
        
        self.sm = numpy.ones(self.b2.shape).astype(bool)
       
        for i in [self.p3m, self.p4m, self.p5m, self.p6m]:
            self.sm = self.CombinM(i, self.sm)
        
        self.eb12 = self.ScaleExtent(self.b12)
        self._del(self.b12)
        b12m = self.eb12 < self.t7_eb12
        
        struct = iterate_structure(generate_binary_structure(2,1), 3)
        snow_mask_dil = binary_dilation(self.sm, struct)
        bz = snow_mask_dil -  self.sm
        bs = bz & b12m
        #print bs
        self.sm = self.CombinM(bs, self.sm, t=False)
        
        self.sp = numpy.ones(self.b2.shape)
        
        for i in [self.p3p, self.p4p, self.p5p, self.p6p]:
            self.sp = self.Posb(i,self.sp)
            
        
        S_p_m = self.sp > self.t7_sp #use the finnal possibility to thershold (0.73 for his) the snow mask
        self.sm = self.CombinM(self.sm, S_p_m)
        
        S_free = ~self.sm
        
        self.cm = self.p1_p2_M & S_free
        
        
        #self.p7p, self.p7m = M_P_func(b2_b4, self.t7_min, self.t7_max)
    def _p8(self):
        
        self.NDVI = (self.b8 - self.b4)/(self.b8 + self.b4)
        p8p, self.p8m = self.M_P_func(self.NDVI, self.t8_min, self.t8_max, nm = 4) 
        self.p8p =  p8p*-1 + 1
        self.vm = ~ self.p8m
        
    
    def _p9(self):
        
        b8_b3 = self.b8/self.b3
        p9p, self.p9m = self.M_P_func(b8_b3, self.t9_min, self.t9_max, nm=4)
        self.p9p = 1. - p9p
        
    def _p10(self):
        
        b2_b11 = self.b2/self.eb11
        
        b2_FT = clip(-0.4*b2_b11 + 0.46, 0.15, 0.32)        
        
        self.bsm = (self.b2 < b2_FT) & (b2_b11 < 0.6) #bare soil
        #self.tempM = (self.b2_b11 > 0.6) & (self.b2_b11 < 0.9) & (self.b2 < b2_FT)
        self.p10p, p10m = self.M_P_func(b2_b11, self.t10_min, self.t10_max)
        self.p10m = p10m & ~(self.b2 < b2_FT)
        #self.b2_FT = b2_FT
        #self.b2_b11 = b2_b11
     
    # water bodies
    def _p11(self):
        
        self.eb8a = self.ScaleExtent(self.b8a)
        d_b2_b4 = self.b2-self.b4
        self.wm = (d_b2_b4 > self.t11_1) & (self.eb8a < self.b4) & (self.b2 < self.t11_2)
        
        self.p11m = ~self.wm 
        self.p11p = 1.
    
    def _p12(self):
        
        b2_b11 = self.b2/self.eb11
        B12_FT = clip(b2_b11*0.1 - 0.09, self.t12_min, self.t12_max)
        self.wm2 = where((b2_b11 > 4.) & (self.eb12 < B12_FT) & (self.eb8a < self.b4) & (self.b2 < 0.2), True, False)
    
        R15_AMB = (b2_b11 < 4.) & (b2_b11 >= 2.) & (self.eb12 < B12_FT) & (self.eb8a < self.b4) & (self.b2 < 0.2)

        #if(R15_AMB.size > 0):
         #   a = -1 / (T22_R_B02_B11 - T21_R_B02_B11)
         #   b = -T21_R_B02_B11 * a + 1
         #   self.p12p = a * R_B02_B11 + b        

        self.wm = R15_AMB | self.wm | self.wm2
        self.p12m = ~self.wm
        self.p12p = 1.
    
    # sand and stone in desert
    def _p13(self):
        
        b8a_b11 = self.eb8a/self.eb11
        b2_FT = clip(b8a_b11*(-0.25)+0.475, 0.16, 0.35)
        self.Dssm = (self.b2 < b2_FT) & (b8a_b11 < 0.9) & (self.b2 < 0.8)
        self.p13p, p13m = self.M_P_func(b8a_b11, self.t13_min, self.t13_max)
        self.p13m = p13m & ~(self.b2 < b2_FT)
        
        
        self.b2_FT = b2_FT
        self.b8a_b11 = b8a_b11
    
    def _p14(self):
        
        b4_b11 = self.b4/self.eb11
        p14p, self.p14m = self.M_P_func(b4_b11, self.t14_min, self.t14_max, nm = 4)
        self.p14p = 1. - p14p
    
    def runit(self, n):
        if n ==1:
            self._p1()
            print 'finished', self.p1p, self.p1m
        

        
import multiprocessing as mp
pool = mp.Pool()

def _p1(test):
    return test._p1()
test.z = pool.apply_async(_p1,args=(test(),)).get()

    def Get_cm_p(self):
        
        self.Getimg(self.bands, self.fhead, self.bounds)
        
        fs = ['self._p%s()'%i for i in range(7)[1:]]
        pool = mp.ProcessingPool(10)
        print pool
        pool.map(self.runit, range(10))
        pool.close()
        pool.join()
        
        print 'finished'
        self._p7()
        
        fs = ['self._p%s()'%i for i in range(7,15)]
        pool = mp.ProcessingPool(10)
        pool.map(self.runit, fs)
        pool.close()
        pool.join()
        '''
        self._p7()
        print 'reading image...'
        self.Getimg(self.fhead, self.bands, self.bounds)
        print 'process 1...'
        self._p1()
        print 'process 2...'
        self._p2()
        print 'process 3...'
        self._p3()
        print 'process 4...'
        self._p4()
        print 'process 5...'
        self._p5()
        print 'process 6...'
        self._p6()
        print 'process 7...'
        self._p7()
        print 'process 8...'
        self._p8()
        print 'process 9...'
        self._p9()
        print 'process 10...'
        self._p10()
        print 'process 11...'
        self._p11()
        print 'process 12...'
        self._p12()
        print 'process 13...'
        self._p13()
        print 'process 14...'
        self._p14()
        '''
        '''
        for i in [self.p1m, sel.p2m,~self.sm, self.p8m, self.p9m, self.p10m, self.p11m, self.p12m, self.p13m, self.p14m]:
            self.cm = self.CombinM(i, self.cm)
        
        self.cp = self.p1_p2_P
        for i in [self.p1p, self.p2p, self.p8p, self.p9p, self.p10p, self.p11p, self.p12p, self.p13p, self.p14p]:
            self.cp = self.Posb(i,self.cp)
        
        c_ptm = self.cp > self.c_pt
        self.fcm = self.cm & c_ptm
        
        v = [self.p8m, self.p9m, self.p10m, self.p11m, self.p12m, self.p13m, self.p14m, 
                  self.p8p, self.p9p, self.p10p, self.p11p, self.p12p, self.p13p, self.p14p]
        self._del(v)
        '''
    
    
'''
{k.strip(): float(v) for k, v in a}
{'self.c_pt': 0.015,
 'self.t10_max': 0.9,
 'self.t10_min': 0.6,
 'self.t11_1': 0.034,
 'self.t11_2': 0.2,
 'self.t12_max': 0.21,
 'self.t12_min': 0.07,
 'self.t13_max': 1.1,
 'self.t13_min': 0.9,
 'self.t14_max': 6.0,
 'self.t14_min': 3.0,
 'self.t1_max': 0.25,
 'self.t1_min': 0.07,
 'self.t2_max': 0.2,
 'self.t2_min': -0.1,
 'self.t3_max': 0.42,
 'self.t3_min': 0.2,
 'self.t4_max': 0.35,
 'self.t4_min': 0.15,
 'self.t5_max': 0.22,
 'self.t5_min': 0.18,
 'self.t6_max': 0.95,
 'self.t6_min': 0.85,
 'self.t7_eb12': 0.35,
 'self.t7_sp': 0.73,
 'self.t8_max': 0.2,
 'self.t8_min': 0.02,
 'self.t9_max': 1.6,
 'self.t9_min': 1.0}
 '''  
def p1(classification):
    return classification._p1()
def p2(classification):
    return classification._p2()
def p3(classification):
    return classification._p3()
def p4(classification):
    return classification._p4()
def p5(classification):
    return classification._p5()
def p6(classification):
    return classification._p6()
def p7(classification):
    return classification._p7()
def p8(classification):
    return classification._p8()
def p9(classification):
    return classification._p9()
def p10(classification):
    return classification._p10()
def p11(classification):
    return classification._p11()
def p12(classification):
    return classification._p12()
def p13(classification):
    return classification._p13()
def p14(classification):
    return classification._p14()


import multiprocessing as mp
pool = mp.Pool(processes=10)

#classification.Getimg(self.bands, self.fhead, self.bounds)

#pool.apply_async(p1,args=(classification(),))      
classification.p1p, classification.p1m = pool.apply_async(p1,args=(classification(img),)).get()
classification.p2p, classification.p2m = pool.apply_async(p2,args=(classification(img),)).get()
classification.p3p, classification.p3m = pool.apply_async(p3,args=(classification(img),)).get()
classification.p4p, classification.p4m = pool.apply_async(p4,args=(classification(img),)).get()
classification.p5p, classification.p5m = pool.apply_async(p5,args=(classification(img),)).get()
classification.p6p, classification.p6m = pool.apply_async(p6,args=(classification(img),)).get()
pool.close()
pool.join()


#fs = ['self._p%s()'%i for i in range(7,15)]
pool = multiprocessing.Pool(processes=10)
self.p8p, self.p8m, self.vm, self.NDVI = pool.apply_async(p8, args=(classification(img),)).get()
self.p9p, self.p9m = pool.apply_async(p9,args=(classification(img),)).get()
self.p10p, self.p10m, self.bsm = pool.apply_async(p10,args=(classification(img),)).get()
self.p11p, self.p11m, self.wm = pool.apply_async(p11, args=(classification(img),)).get()
self.p12p, self.p12m, self.wm2,self.wm3 = pool.apply_async(p12,args=(classification(img),)).get()
self.p13p, self.p13m, self.Dssm = pool.apply_async(p13, args=(classification(img),)).get()
self.p14p,self.p14m = pool.apply_async(p14, args=(classification(img),)).get()
pool.close()
pool.join()
class test(object):
    z =10
    def p1(self):
        
        return self.z
    def p2(self):
        
        
        z = 10**6
        
    def runit(self):
        pool = ProcessingPool(nodes=20)
        pool.map(self.p1,range(100,103))
        pool.map(self.p2,range(1000000,1000001))
import multiprocessing as mp
pool = mp.Pool()


def p1(test):
    return test.p1()

pool.apply_async(p1,args=(test(),)).get()
"""        
