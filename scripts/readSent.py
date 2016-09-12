import glymur 
import gdal
import xml.etree.ElementTree as ET
import numpy 

'''
    Read in all of the data from sentinel2 L1C files
    
    ---Marc Yin
    24/05/2016

'''

def readfile(fhead,band=xrange(1,13),pr=False):
    '''
    This function read in all of the img data, meteological data,
    and the sun and viewing geometry data, apart from the gml file
    it returns all a dict store the keys and the array of the data.
    
    parameters:
    fhead--filehead like '50SMG20164100'(Military Grid Reference System tile, doy,and '0'(subdirectory))
    band-- default to be all of the bands (1-12 and A8 band)
    pr-- defult to be False, if True then print the jp2 file information, product information and tile information
    
    return:
    a dict of different arrays:
    
    TCWV: total column water vapor
    MSLP: mean sea level presure
    TCO3: total column ozone 
    SAG--Sun_Angles_Grid 
    VIAG--Viewing_Incidence_Angles_Grids
    
    mSz,mSa,mVz,mVa--
        m--mean
        a--azimuth
        z--zenith
    
    -- Marc Yin
    24/05/2016
    '''
    
    imgdata = readimg('data/'+fhead,band,pr)
    if pr:
        print 'tileInfo: ', open('data/%stileInfo.json'%fhead,'rb').read(),'\n',
        'prodInfo: ', open('data/%sproductInfo.json'%fhead,'rb').read()
    
    meteo = readmeteo('data/%sauxiliaryECMWFT'%fhead)
    xmlfile = readxml('data/%smetadata.xml'%fhead)
    
    imgdata.update(xmlfile)
    imgdata.update(meteo)
    
    return imgdata
    
def readimg(fhead,band,pr = False):
    filenames = []
    for i in band:
        if i ==8:
            e = fhead+'B0'+'%s'%i +'.jp2'
            filenames.append(e)
            e1 = fhead+'B'+'%sA'%i +'.jp2'
            filenames.append(e1)
        elif i<10:
            e = fhead+'B0'+'%s'%i +'.jp2'
            filenames.append(e)
        else:
            filename = fhead+'B'+'%s'%i +'.jp2'
            filenames.append(filename)
    #print filenames

    imgdata = {k.split('.')[0][-3:]: [] for k in filenames}

    for i in filenames:
        jp2file = glymur.data.nemo()
        #os.chdir('data')
        jp2 = glymur.Jp2k(i)
        fullres = jp2[:]
        if pr:
            print jp2,'\n',fullres.shape
        
        imgdata[i.split('.')[0][-3:]] = fullres/10000. 
    
    return imgdata

def readmeteo(filename):
    '''
    read in meteological data
    TCWV: total column water vapor
    MSLP: mean sea level presure
    TCO3: total column ozone 
    '''
    g = gdal.Open(filename)
    if g is None:
        print 'No meteo data!!!'
    data = g.ReadAsArray()
    return {'TCWV':data[0],'MSLP':data[1],'TCO3':data[2]}




def readxml(filename):
    '''
    This function is only used for the Sentinel 2 L1C metadata.xml file
    the major function of this module is to get the sun zenith angle and viewing angle
    grid and mean value are provided with the coloum and raw having a step value of 5k m
    for the grid it has 13 band (0-12) and each band have 12 detector_id
    
    in:
    filename: with the right path!!
    out:
    a dict: use the dict.keys() to check the file key names
    the abbrivation: SAG--Sun_Angles_Grid; VIAG--Viewing_Incidence_Angles_Grids; m--mean; A--Azimuth; Z--Zenith
    ---Marc Yin
    23/05/2016
    '''
   
    
    tree = ET.parse(filename)
    root = tree.getroot()
    #Sun_Angles_Grid
    SAG_A=[]
    SAG_Z=[]
    mSz = []
    mSa = []
    #Viewing_Incidence_Angles_Grids
    VIAG_A = []
    VIAG_Z = []
    mVz = []
    mVa = []
    for child in root:
        for j in child:

            for k in j.findall('Sun_Angles_Grid'):
                for l in k.findall('Zenith'):
                    for m in l.findall('Values_List'):
                        for x in m.findall('VALUES'):
                            SAG_A.append(x.text)
                for n in k.findall('Azimuth'):
                    for o in n.findall('Values_List'):
                        for p in o.findall('VALUES'):
                            SAG_Z.append(p.text) 

            for msa in j.findall('Mean_Sun_Angle'):
                mSz.append(msa.find('ZENITH_ANGLE').text)
                mSa.append(msa.find('AZIMUTH_ANGLE').text)

            #for viag in j.iter('Viewing_Incidence_Angles_Grids'):
            for k in j.findall('Viewing_Incidence_Angles_Grids'):
                for l in k.findall('Zenith'):
                    for m in l.findall('Values_List'):
                        for x in m.findall('VALUES'):
                            VIAG_A.append(x.text)
                for n in k.findall('Azimuth'):
                    for o in n.findall('Values_List'):
                        for p in o.findall('VALUES'):
                            VIAG_Z.append(p.text)

            for mvia in j.findall('Mean_Viewing_Incidence_Angle_List'):
                for i in mvia.findall('Mean_Viewing_Incidence_Angle'):
                    mVz.append(i.find('ZENITH_ANGLE').text)
                    mVa.append(i.find('AZIMUTH_ANGLE').text)




    SAG_A = [(i.split(' ')) for i in SAG_A]
    SAG_A = numpy.array(SAG_A).astype(float)
    SAG_Z = [(i.split(' ')) for i in SAG_Z]
    SAG_Z = numpy.array(SAG_Z).astype(float)
    mSa = numpy.array(mSa).astype(float)
    mSz = numpy.array(mSz).astype(float)
    shape = (len([i.split(' ') for i in VIAG_A])/23,23,23)        
    VIAG_A = numpy.array([i.split(' ') for i in VIAG_A]).reshape(shape).astype(float)
    shape = (len([i.split(' ') for i in VIAG_Z])/23,23,23)
    VIAG_Z = numpy.array([i.split(' ') for i in VIAG_Z]).reshape(shape).astype(float)
    mVa = numpy.array(mVa).astype(float)
    mVz = numpy.array(mVz).astype(float)
    
    return {'SAG_A':SAG_A,'SAG_Z':SAG_Z,'mSa':mSa,'mSz':mSz,
            'VIAG_A':VIAG_A,'VIAG_Z':VIAG_Z,'mVa':mVa,'mVz':mVz}
    
    
   
