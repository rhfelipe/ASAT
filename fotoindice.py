from PyQt6.QtCore import QThread, Qt
from PyQt6.QtCore import QThread, pyqtSignal, pyqtSlot, QRegularExpression
import sys
import rasterio as rst
import rasterio.sample
import rasterio.vrt
import rasterio._features
from PIL import Image, ImageDraw, ImageFont
import os
import numpy as np
import math as mt
from rasterio.transform import Affine
import csv
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
import matplotlib.ticker as ticker
import time
from pyproj import Transformer
import geopandas as gpd
from fiona.drvsupport import supported_drivers
from io import BytesIO
from reportlab.lib.pagesizes import letter, A3
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
import tempfile
from matplotlib.backends.backend_pdf import FigureCanvasPdf
from reportlab.lib.utils import ImageReader
import utm
import ee

diretorio_atual = os.path.abspath(os.path.dirname(sys.executable))
#diretorio_atual = os.path.abspath(os.path.dirname(__file__))


ee_account='conta google earth engine'
credentials = ee.ServiceAccountCredentials(ee_account, diretorio_atual + '/resources/private-key.json')
ee.Initialize(credentials)

class Thread_1(QThread):
    # Sinais para comunicação entre as threads
    finished = pyqtSignal(int)
    result_ready = pyqtSignal(list)
    progress_changed = pyqtSignal(int, str, bool, int)
    etapa_changed = pyqtSignal(int)

    def __init__(self, tifs, dir, pixel_size, fuso):
        super().__init__()
        self.tifs = tifs
        self.dir = dir
        self.pixel_size = int(pixel_size)/100
        self.fuso = fuso


    def run(self):
        # Lógica da thread
        result = self.imagens_coords()
        # Emitir o sinal com o resultado
        if result:
           self.result_ready.emit(result)

        # Emitir o sinal de finalização
    def imagens_coords(self):
        _, arq = os.path.split(self.dir)
        _, typ = arq.split('.')
        if typ == 'csv':
                xn = []
                yn = []
                file = open(self.dir,'r')
                reader = csv.reader(file)
                tam_reader = len(list(reader))
                imagens_coords = {}
                i = 0
                file.seek(0)
                reader = csv.reader(file)
                for row in reader:
                    infos = row[0].split(' ')
                    if len(infos)>=7:
                        xn.append(infos[1])
                        yn.append(infos[2])
                        Z = infos[3].replace(',','.')
                        Z = float(Z.replace(" m",""))

                        pitch = mt.radians(float(infos[9].replace(',','.')))
                        roll = mt.radians(float(infos[10].replace(',','.')))

                        imagens_coords[infos[0]]=[infos[1],infos[2], infos[6], infos[3]]
                        self.progress_changed.emit((i*100/tam_reader)/2,"",False,1)
                        i=i+1
                    else:
                        self.progress_changed.emit(150,"",True,1)
                        #sys.exit()
                        break
    
                if len(self.tifs)!=len(imagens_coords):
                    tempo = {}
                    xn = []
                    yn = []
                    tam_tifs = len(self.tifs)
                    i = 0
                    for tf in self.tifs:
                        _,nm = os.path.split(tf)
                        try:
                           vet = imagens_coords[nm]
                           tempo[nm]=vet
                           xn.append(vet[0])
                           yn.append(vet[1])
                        except:
                            pass
                        self.progress_changed.emit(50+(i*100/tam_tifs)/2,"",False,1)
                        i=i+1
                    imagens_coords = {}
                    imagens_coords = tempo
    
    
                if len(imagens_coords)>0:
                    sizes = [float(max(xn)),float(min(xn)),float(max(yn)),float(min(yn))]
                    dif_x = int(sizes[0] - sizes[1])
                    dif_y = int(sizes[2] - sizes[3])
    
                    img_temp = Image.open(self.tifs[0])
                    og_w = img_temp.width
                    og_h = img_temp.height
                    img_temp.close()
    
                    if dif_x >= dif_y:
                        scale = dif_x/406.04694 
                        ex_x = int((dif_x/scale)/(self.pixel_size*10)) #multiplicando por 10 o gsd para que seja o tamanho da foto no A3
                        fator = og_w/ex_x
                        ex_y = int(og_h/fator)
                        dist_y=(279.48*scale)
                        acrescimo = mt.fabs(dif_y-dist_y)/2
                        bds = [float(max(xn))+2000, float(min(xn))-2000, float(max(yn))+acrescimo, float(min(yn))-acrescimo]
                        ind = 1
                    else:
                        scale = dif_y/406.04694 
                        ex_y = int((dif_y/scale)/(self.pixel_size*10))
                        fator = og_h/ex_y
                        ex_x = int(og_w/fator)
                        dist_x=(279.48*scale)
                        acrescimo = mt.fabs(dif_x-dist_x)/2
                        bds = [float(max(xn))+acrescimo, float(min(xn))-acrescimo, float(max(yn))+2000, float(min(yn))-2000]
                        ind = 2
                    size = [ex_x*2,ex_y*2]
                    self.etapa_changed.emit(1)
                    #self.finished.emit(1)
                    return [imagens_coords, size, bds, scale, ind]
                else:
                    self.progress_changed.emit(152,"",True,1)
    
                file.close()
        elif typ == 'txt':
                def images():
                         #Imagens
                         imagem_2 = ee.Image('CGIAR/SRTM90_V4')
                         elevacao = imagem_2.select('elevation')
                         slope = ee.Terrain.slope(elevacao)
                         topo = ee.Image.cat(elevacao, slope)
                         topocol = ee.ImageCollection([topo])

                         return topocol
                
                def rasterExtraction(image):
                              try:
                                 feature = image.sample(region = point, geometries = False)
                                 return feature
                              except:
                                  return False


                with open(self.dir, 'r') as arquivo:
                    # Ler todas as linhas do arquivo
                    linhas = arquivo.readlines()
                    linhas = linhas[1:]
                    i = 0
                    xn = []
                    yn = []
                    imagens_coords = {}
                    tam_reader = len(linhas)
                    _, tft = os.path.split(self.tifs[0])
                    tptft = tft.split('.')[1]

                    topocol = images()
                    for linha in linhas:
                        infos = linha.split("	")
                        if 'X' in infos[6] or not infos[6]  or 'X' in infos[7] or not infos[7] or 'Longitude' in infos[6]:
                            continue
                        else:
                            if tptft == 'tif':
                                foto = infos[0].replace(".IIQ",".tif")
                            else:
                                foto = infos[0]
                            long = float(infos[6].replace(",","."))
                            lat = float(infos[7].replace(",","."))
                            yaw = float(infos[11].replace(",","."))

                            point = ee.Geometry.Point(long, lat)

                            results = topocol.filterBounds(point).map(rasterExtraction).flatten()
                            if results:
                                      resultados = results.getInfo()

                                      if len(resultados['features'])==0 and 'elevation' in locals():
                                          elevation = elevation
                                      elif len(resultados['features'])==0 and 'elevation' not in locals():
                                          elevation=0
                                      else:
                                          elevation = resultados['features'][0]['properties']['elevation']

                                      
                            else:
                                elevation = 0





                            Z = infos[8].replace(',','.')
                            Z = float(Z.replace(" m",""))-elevation

                            pitch = mt.radians(float(infos[9].replace(',','.')))
                            roll = mt.radians(float(infos[10].replace(',','.')))

                            

                           


                            en = list(utm.from_latlon(lat, long, int(self.fuso)))
                            xn.append(en[0])
                            yn.append(en[1])
                            imagens_coords[foto]=[en[0],en[1], yaw, Z, pitch, roll]
                            self.progress_changed.emit(int((i*100/tam_reader)/2),"",False,1)
                            i=i+1

                    if len(self.tifs)!=len(imagens_coords):
                        tempo = {}
                        xn = []
                        yn = []
                        tam_tifs = len(self.tifs)
                        i = 0
                        for tf in self.tifs:
                            _,nm = os.path.split(tf)
                            try:
                               vet = imagens_coords[nm]
                               tempo[nm]=vet
                               xn.append(vet[0])
                               yn.append(vet[1])
                            except:
                                pass
                            self.progress_changed.emit(int(50+(i*100/tam_tifs)/2),"",False,1)
                            i=i+1
                        imagens_coords = {}
                        imagens_coords = tempo

                    if len(imagens_coords)>0:
                        sizes = [float(max(xn)),float(min(xn)),float(max(yn)),float(min(yn))]
                        dif_x = int(sizes[0] - sizes[1])
                        dif_y = int(sizes[2] - sizes[3])
        
                        img_temp = Image.open(self.tifs[0])
                        og_w = img_temp.width
                        og_h = img_temp.height
                        img_temp.close()
        
                        if dif_x >= dif_y:
                            scale = dif_x/406.04694 # pensando numa folha A3 deitada + margem
                            ex_x = int((dif_x/scale)/(self.pixel_size*10)) #multiplicando por 10 o gsd para que seja o tamanho da foto no A3
                            fator = og_w/ex_x
                            ex_y = int(og_h/fator)
                            dist_y=(279.48*scale)
                            acrescimo = mt.fabs(dif_y-dist_y)/2
                            bds = [float(max(xn))+3000, float(min(xn))-3000, float(max(yn))+acrescimo, float(min(yn))-acrescimo]
                            ind = 1
                        else:
                            scale = dif_y/406.04694 #pensando numa folha A3 em pé + margem
                            ex_y = int((dif_y/scale)/(self.pixel_size*10))
                            fator = og_h/ex_y
                            ex_x = int(og_w/fator)
                            dist_x=(279.48*scale)
                            acrescimo = mt.fabs(dif_x-dist_x)/2
                            bds = [float(max(xn))+acrescimo, float(min(xn))-acrescimo, float(max(yn))+3000, float(min(yn))-3000]
                            ind = 2
                        size = [ex_x*2,ex_y*2]
                        self.etapa_changed.emit(1)
                        #self.finished.emit(1)
                        return [imagens_coords, size, bds, scale, ind]
                    else:
                        self.progress_changed.emit(152,"",True,1)



#foi removido o crs=self.epsg pois quando convertido ao exe estava dando algum problema
class Thread_2(QThread):
    # Sinais para comunicação entre as threads
    finished = pyqtSignal()
    result_ready = pyqtSignal(list, int)
    progress_changed = pyqtSignal(int, str, bool, int)
    etapa_changed = pyqtSignal(int)

    def __init__(self, tifs, size, fuso, gsd, imagens_coords, id):
        super().__init__()
        self.tifs = tifs
        self.size = size
        self.epsg = 'EPSG:' + str(31960+int(fuso))
        self.pixel_size = gsd/100
        self.imagens_coords = imagens_coords
        img_temp = Image.open(tifs[0])
        self.og_w = img_temp.width
        self.og_h = img_temp.height
        self.id = id
        self.fuso = int(fuso)

        img_temp.close()


    def run(self):
        # Lógica da thread
        result = self.mem_file()

        # Emitir o sinal com o resultado
        self.result_ready.emit(result, 1)

        # Emitir o sinal de finalização
        self.finished.emit()

    def coords_cantos(self, coords, altura, largura, modelo, fuso):
         gsd = self.pixel_size
         altura_o = self.og_h/2
         largura_o = self.og_w/2
         ang = float(coords[2])
         z = coords[3]
         if largura_o < 1000:
             if modelo == 'iXM-RS100F' or modelo=='iXU-RS1000':
                 gsd = (z/51.7383)*0.00534
                 gsd = gsd*18.1375
             else:
                 gsd = (z/45.9310)*0.00537
                 gsd = gsd*16.1375
         ly = altura_o*gsd
         lx = largura_o*gsd
         roll = float(coords[5])
         pitch = float(coords[4])
         dy = mt.sin(mt.radians(pitch))*ly
         dx= mt.sin(mt.radians(roll))*lx

         x = float(coords[0])
         y = float(coords[1])
         ang = float(coords[2])
         dist = mt.sqrt(((altura_o)*gsd)**2 + ((largura_o)*gsd)**2)
         pgcp1 = list(utm.to_latlon(x + dist*mt.sin((315+ang)*mt.pi/180)-dx, y + dist*mt.cos((315+ang)*mt.pi/180)-dy, fuso, northern=False))
         pgcp2 = list(utm.to_latlon(x + dist*mt.sin((ang+45)*mt.pi/180)-dx, y + dist*mt.cos((45+ang)*mt.pi/180)-dy, fuso, northern=False))
         pgcp3 = list(utm.to_latlon(x, y, fuso, northern=False))
         pgcp4 = list(utm.to_latlon(x + dist*mt.sin((ang+135)*mt.pi/180)-dx, y + dist*mt.cos((135+ang)*mt.pi/180)-dy, fuso, northern=False))
         pgcp5 = list(utm.to_latlon(x + dist*mt.sin((ang+225)*mt.pi/180)-dx, y + dist*mt.cos((225+ang)*mt.pi/180)-dy, fuso, northern=False))
         gcp1 = rst.control.GroundControlPoint(0, 0, pgcp1[1], pgcp1[0])
         gcp2 = rst.control.GroundControlPoint(0, largura, pgcp2[1], pgcp2[0])
         gcp3 = rst.control.GroundControlPoint(altura/2, largura/2, pgcp3[1], pgcp3[0])
         gcp4 = rst.control.GroundControlPoint(altura, largura, pgcp4[1], pgcp4[0])
         gcp5 = rst.control.GroundControlPoint(altura, 0, pgcp5[1], pgcp5[0])
         return [gcp1, gcp2, gcp3, gcp4, gcp5]

    def mem_file(self):

        geotiffs = []
        tam_tifs = len(self.tifs)
        i=0
        size = self.size
        imagens_coords = self.imagens_coords
        for tif in self.tifs:
            img = Image.open(tif)

            #if img.width > 1500:
             #  img.thumbnail((size[0],size[1]))
            exif = img.getexif()
            img_np = np.array(img)
            modelo = exif.get(272)

            _, names = os.path.split(tif)

            try:
               xy = imagens_coords[names]
               altura, largura = img_np.shape[:2]
               num_bandas = img_np.shape[2]

               #transform = Affine.translation(float(xy[0]) - (size[0]*pixel_size/2), float(xy[1]) + (size[1]*pixel_size/2))
               #transform.rotation(float(xy[2]))
               gcps = self.coords_cantos(xy, altura, largura, modelo, self.fuso)
               transform = rst.transform.from_gcps(gcps)
               img_np = np.transpose(img_np, (2, 0, 1))


               img.close()

               memfile = rst.io.MemoryFile()
               crs=rst.CRS.from_wkt('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]')
               with memfile.open(driver='GTiff',width=largura, height=altura, count=num_bandas, dtype=img_np.dtype,transform=transform,crs=crs) as dst:
                   dst.write(img_np)
               geotiffs.append(memfile)
               self.progress_changed.emit(int((i*100/tam_tifs)),names,False,2)
            except Exception as e:
               with open(diretorio_atual+'/resources/log.txt', 'w') as f:
                    f.write(str(e))
               self.progress_changed.emit(int((i*100/tam_tifs)),names,True,2)
            i=i+1
        return geotiffs

class Thread_3(QThread):
    # Sinais para comunicação entre as threads
    finished = pyqtSignal()
    result_ready = pyqtSignal(list, int)
    progress_changed = pyqtSignal(int, str, bool, int)
    etapa_changed = pyqtSignal(int)

    def __init__(self, tifs, size, fuso, gsd, imagens_coords, id):
        super().__init__()
        self.tifs = tifs
        self.size = size
        self.epsg = 'EPSG:' + str(31960+int(fuso))
        self.pixel_size = gsd/100
        self.imagens_coords = imagens_coords
        img_temp = Image.open(tifs[0])
        self.og_w = img_temp.width
        self.og_h = img_temp.height
        self.id = id
        self.fuso = int(fuso)

        img_temp.close()


    def run(self):
        # Lógica da thread
        result = self.mem_file()

        # Emitir o sinal com o resultado
        self.result_ready.emit(result, 2)

        # Emitir o sinal de finalização
        self.finished.emit()

    def coords_cantos(self, coords, altura, largura, modelo, fuso):
         gsd = self.pixel_size
         altura_o = self.og_h/2
         largura_o = self.og_w/2
         ang = float(coords[2])
         z = coords[3]
         if largura_o < 1000:
             if modelo == 'iXM-RS100F' or modelo=='iXU-RS1000':
                 gsd = (z/51.738)*0.00534
                 gsd = gsd*18.1375
             else:
                 gsd = (z/45)*0.00537
                 gsd = gsd*16.1375
         ly = altura_o*gsd
         lx = largura_o*gsd
         roll = float(coords[5])
         pitch = float(coords[4])
         dy = mt.sin(mt.radians(pitch))*ly
         dx= mt.sin(mt.radians(roll))*lx
         x = float(coords[0])
         y = float(coords[1])
         ang = float(coords[2])
         dist = mt.sqrt(((altura_o)*gsd)**2 + ((largura_o)*gsd)**2)
         pgcp1 = list(utm.to_latlon(x + dist*mt.sin((315+ang)*mt.pi/180)-dx, y + dist*mt.cos((315+ang)*mt.pi/180)-dy, fuso, northern=False))
         pgcp2 = list(utm.to_latlon(x + dist*mt.sin((ang+45)*mt.pi/180)-dx, y + dist*mt.cos((45+ang)*mt.pi/180)-dy, fuso, northern=False))
         pgcp3 = list(utm.to_latlon(x, y, fuso, northern=False))
         pgcp4 = list(utm.to_latlon(x + dist*mt.sin((ang+135)*mt.pi/180)-dx, y + dist*mt.cos((135+ang)*mt.pi/180)-dy, fuso, northern=False))
         pgcp5 = list(utm.to_latlon(x + dist*mt.sin((ang+225)*mt.pi/180)-dx, y + dist*mt.cos((225+ang)*mt.pi/180)-dy, fuso, northern=False))
         gcp1 = rst.control.GroundControlPoint(0, 0, pgcp1[1], pgcp1[0])
         gcp2 = rst.control.GroundControlPoint(0, largura, pgcp2[1], pgcp2[0])
         gcp3 = rst.control.GroundControlPoint(altura/2, largura/2, pgcp3[1], pgcp3[0])
         gcp4 = rst.control.GroundControlPoint(altura, largura, pgcp4[1], pgcp4[0])
         gcp5 = rst.control.GroundControlPoint(altura, 0, pgcp5[1], pgcp5[0])
         return [gcp1, gcp2, gcp3, gcp4, gcp5]

    def mem_file(self):

        geotiffs = []
        tam_tifs = len(self.tifs)
        i=0
        size = self.size
        imagens_coords = self.imagens_coords
        for tif in self.tifs:
            img = Image.open(tif)

            #img.thumbnail((size[0],size[1]))

            exif = img.getexif()
            img_np = np.array(img)
            modelo = exif.get(272)

            img_np = np.array(img)

            _, names = os.path.split(tif)

            try:
               xy = imagens_coords[names]
               altura, largura = img_np.shape[:2]
               num_bandas = img_np.shape[2]

               #transform = Affine.translation(float(xy[0]) - (size[0]*pixel_size/2), float(xy[1]) + (size[1]*pixel_size/2))
               #transform.rotation(float(xy[2]))
               gcps = self.coords_cantos(xy, altura, largura, modelo, self.fuso)
               transform = rst.transform.from_gcps(gcps)
               img_np = np.transpose(img_np, (2, 0, 1))


               img.close()

               memfile = rst.io.MemoryFile()
               crs=rst.CRS.from_wkt('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]')
               with memfile.open(driver='GTiff',width=largura, height=altura, count=num_bandas, dtype=img_np.dtype,transform=transform,crs=crs) as dst:
                   dst.write(img_np)
               geotiffs.append(memfile)
               self.progress_changed.emit(int((i*100/tam_tifs)),names,False,3)
            except Exception as e:
               self.progress_changed.emit(int((i*100/tam_tifs)),names,True,3)
            i=i+1
        return geotiffs


class Thread_4(QThread):
    # Sinais para comunicação entre as threads
    finished = pyqtSignal()
    result_ready = pyqtSignal(list, int)
    progress_changed = pyqtSignal(int, str, bool, int)
    etapa_changed = pyqtSignal(int)

    def __init__(self, tifs, size, fuso, gsd, imagens_coords, id):
        super().__init__()
        self.tifs = tifs
        self.size = size
        self.epsg = 'EPSG:' + str(31960+int(fuso))
        self.pixel_size = gsd/100
        self.imagens_coords = imagens_coords
        img_temp = Image.open(tifs[0])
        self.og_w = img_temp.width
        self.og_h = img_temp.height
        self.id = id
        self.fuso = int(fuso)

        img_temp.close()


    def run(self):
        # Lógica da thread
        result = self.mem_file()

        # Emitir o sinal com o resultado
        self.result_ready.emit(result, 3)

        # Emitir o sinal de finalização
        self.finished.emit()

    def coords_cantos(self, coords, altura, largura, modelo, fuso):
         gsd = self.pixel_size
         altura_o = self.og_h/2
         largura_o = self.og_w/2
         ang = float(coords[2])
         z = coords[3]
         if largura_o < 1000:
             if modelo == 'iXM-RS100F' or modelo=='iXU-RS1000':
                 gsd = (z/51.738)*0.00534
                 gsd = gsd*18.1375
             else:
                 gsd = (z/45)*0.00537
                 gsd = gsd*16.1375
         ly = altura_o*gsd
         lx = largura_o*gsd
         roll = float(coords[5])
         pitch = float(coords[4])
         dy = mt.sin(mt.radians(pitch))*ly
         dx= mt.sin(mt.radians(roll))*lx
         x = float(coords[0])
         y = float(coords[1])
         ang = float(coords[2])
         dist = mt.sqrt(((altura_o)*gsd)**2 + ((largura_o)*gsd)**2)
         pgcp1 = list(utm.to_latlon(x + dist*mt.sin((315+ang)*mt.pi/180)-dx, y + dist*mt.cos((315+ang)*mt.pi/180)-dy, fuso, northern=False))
         pgcp2 = list(utm.to_latlon(x + dist*mt.sin((ang+45)*mt.pi/180)-dx, y + dist*mt.cos((45+ang)*mt.pi/180)-dy, fuso, northern=False))
         pgcp3 = list(utm.to_latlon(x, y, fuso, northern=False))
         pgcp4 = list(utm.to_latlon(x + dist*mt.sin((ang+135)*mt.pi/180)-dx, y + dist*mt.cos((135+ang)*mt.pi/180)-dy, fuso, northern=False))
         pgcp5 = list(utm.to_latlon(x + dist*mt.sin((ang+225)*mt.pi/180)-dx, y + dist*mt.cos((225+ang)*mt.pi/180)-dy, fuso, northern=False))
         gcp1 = rst.control.GroundControlPoint(0, 0, pgcp1[1], pgcp1[0])
         gcp2 = rst.control.GroundControlPoint(0, largura, pgcp2[1], pgcp2[0])
         gcp3 = rst.control.GroundControlPoint(altura/2, largura/2, pgcp3[1], pgcp3[0])
         gcp4 = rst.control.GroundControlPoint(altura, largura, pgcp4[1], pgcp4[0])
         gcp5 = rst.control.GroundControlPoint(altura, 0, pgcp5[1], pgcp5[0])
         return [gcp1, gcp2, gcp3, gcp4, gcp5]

    def mem_file(self):

        geotiffs = []
        tam_tifs = len(self.tifs)
        i=0
        size = self.size
        imagens_coords = self.imagens_coords
        for tif in self.tifs:
            img = Image.open(tif)

            #if img.width > 1500:
               #img.thumbnail((size[0],size[1]))

            exif = img.getexif()
            img_np = np.array(img)
            modelo = exif.get(272)

            img_np = np.array(img)

            _, names = os.path.split(tif)

            try:
               xy = imagens_coords[names]
               altura, largura = img_np.shape[:2]
               num_bandas = img_np.shape[2]

               #transform = Affine.translation(float(xy[0]) - (size[0]*pixel_size/2), float(xy[1]) + (size[1]*pixel_size/2))
               #transform.rotation(float(xy[2]))
               gcps = self.coords_cantos(xy, altura, largura, modelo, self.fuso)
               transform = rst.transform.from_gcps(gcps)
               img_np = np.transpose(img_np, (2, 0, 1))


               img.close()

               memfile = rst.io.MemoryFile()
               crs=rst.CRS.from_wkt('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]')
               with memfile.open(driver='GTiff',width=largura, height=altura, count=num_bandas, dtype=img_np.dtype,transform=transform,crs=crs) as dst:
                   dst.write(img_np)
               geotiffs.append(memfile)
               self.progress_changed.emit(int((i*100/tam_tifs)),names,False,4)
            except Exception as e:
               self.progress_changed.emit(int((i*100/tam_tifs)),names,True,4)
            i=i+1
        return geotiffs

class Thread_5(QThread):
    # Sinais para comunicação entre as threads
    finished = pyqtSignal()
    result_ready = pyqtSignal()
    progress_changed = pyqtSignal(int, str, bool, int)
    etapa_changed = pyqtSignal(int)

    def __init__(self, geotiffs, kml, pasta, projeto, data, gsd, fuso, escala, epsg, bds, ind, dir, tifs, coords):
        super().__init__()
        self.geotiffs = geotiffs
        self.kml = kml
        self.projeto = projeto
        self.data = data
        self.gsd = gsd
        self.fuso = fuso
        self.escala = escala
        self.epsg = epsg
        self.bds = bds
        self.pasta = pasta
        self.ind = ind
        self.dir = dir
        self.tifs = tifs
        self.coords = coords


    def run(self):
        # Lógica da thread
        result = self.fig_fotoindice()

        # Emitir o sinal com o resultado
        self.result_ready.emit()

        # Emitir o sinal de finalização
        self.finished.emit()

    def cria_kml(self, dic_kml,epsg):
        supported_drivers['KML'] = 'rw'
        data = gpd.read_file(dic_kml,driver='kml')

        #data = data.to_crs(epsg.split(':')[1])


        return data

    def fig_fotoindice(self):
        bds = self.bds
        fig, ax = plt.subplots()
        if self.ind == 1:
            width, height = A3
        else:
            height, width = A3
        #fig, ax= plt.subplots(figsize=(r1, r2))

        # Carregue a imagem do gráfico a partir do objeto BytesIO

        width_i = width-32#16.53 * 1000 # Largura do gráfico no documento PDF
        height_i =  height-32#11.69 * 1000 # Altura do gráfico no documento PDF

        fig.set_size_inches((width_i/inch,height_i/inch))  # Defina o tamanho da figura em polegadas (largura, altura)
        tam_geotifs = len(self.geotiffs)
        i=0
        #plt.figure(figsize=(9.7, 6.3))  # Largura: 8 polegadas, Altura: 6 polegadas
        if self.kml:
            qkml = self.cria_kml(self.kml, self.epsg)
            qkml.plot(ax=ax, facecolor='none', edgecolor='red', linewidth=0.4)
        for geotiff_path in self.geotiffs:
            dataset = geotiff_path.open()
            extent = dataset.bounds  # Obtém os limites da imagem
            transform = dataset.transform


            image = dataset.read()


            #image_data = transparent_img(image,imagens_coords,geotiff_path)
            #image = np.transpose(image, (2,0,1))
            image = np.transpose(image, (1,2,0))

            # Rotacionar a imagem
            angle = np.degrees(mt.atan2(transform[1], transform[0]))
            rotated_image = rotate(image, angle, reshape=True)

            rotated_image = Image.fromarray(rotated_image)
            imagem_rotacionada = rotated_image.convert('RGBA')
            imagem_rotacionada = np.array(imagem_rotacionada)

            pixels_pretos = np.all(imagem_rotacionada[:,:,:3] == [0, 0, 0], axis=-1)
            imagem_rotacionada[pixels_pretos] = [0, 0, 0, 0]

            self.progress_changed.emit(int((i*100/tam_geotifs)),"",False, 5)
            i = i+1



            ax.imshow(imagem_rotacionada,extent=[extent.left,extent.right,extent.bottom,extent.top],alpha=imagem_rotacionada[:,:,3]/255)

        formatter = ticker.StrMethodFormatter('{x:.2f}') 
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        plt.tick_params(axis='both', direction='in')
        plt.gca().xaxis.set_major_formatter(formatter)
        plt.gca().yaxis.set_major_formatter(formatter)
        try:
            pts_max = list(utm.to_latlon(bds[0], bds[2], int(self.fuso) , northern=False))
            pts_min = list(utm.to_latlon(bds[1], bds[3], int(self.fuso), northern=False))
            ax.set_xlim(pts_min[1], pts_max[1])
            plt.xticks(fontsize=8)
            ax.set_ylim(pts_min[0], pts_max[0])
            plt.yticks(fontsize=8,rotation=90,verticalalignment='center')
            #plt.grid(linewidth=0.1,color='black')
            #plt.show()
            buffer = BytesIO()
            plt.savefig(buffer,format='png', dpi=1000,bbox_inches='tight', pad_inches=0.01)
            buffer.seek(0)
    
            # ADICIONANDO O PLOT NO PDF
    
            # Abrir a imagem com PIL
            imagem_pil = Image.open(buffer)
            # Agora você pode trabalhar com a imagem usando a biblioteca PIL, por exemplo:
            larg, alt = imagem_pil.size
    
            if larg >= alt:
                fator = larg/width_i
                w = larg/fator
                h = alt/fator
            else:
                fator = alt/height_i
                h = alt/fator
                w = larg/fator
    
    
             #criando infos adicionais
            img_t = Image.open(self.tifs[0])
            exif = img_t.getexif()
            camera = exif.get(272)
            total_fotos = str(i+1)
            exposicao = '1/1000s'
            img_t.close()
    
            lat = []
            long = []
            z = []
            for key in self.coords:
                lis = self.coords[key]
                lat.append(lis[0])
                long.append(lis[1])
                z.append(lis[3])
            mean_l = list(utm.to_latlon(np.mean(lat), np.mean(long), int(self.fuso), northern=False))
            lat_med = np.format_float_positional(mean_l[0],precision=4)
            long_med = np.format_float_positional(mean_l[1],precision=4)
            z_med = np.format_float_positional(np.mean(z),precision=3)
    
    
             #ABRIR IMAGEM DE INFO DO RODAPÉ
            infos = Image.open(self.dir+'/resources/infos_2.png')
            desenho = ImageDraw.Draw(infos)
            fonte = ImageFont.truetype("arial.ttf", 40)
    
            larg_info, alt_info = infos.size
    
            desenho.text((1880,50), self.projeto ,fill=(0, 0, 0),font=fonte)
            desenho.text((1790,420), self.data ,fill=(0, 0, 0),font=fonte)
            desenho.text((1950,320), str(self.gsd)+'cm',fill=(0, 0, 0),font=fonte)
            desenho.text((1540,320), 'WGS84/' + str(self.fuso) + 'S' ,fill=(0, 0, 0),font=fonte)
            desenho.text((1490,420), '1/'+ str(int(self.escala*1000)) ,fill=(0, 0, 0),font=fonte)
            desenho.text((815,205),total_fotos ,fill=(0, 0, 0),font=fonte)
            desenho.text((1120,205),camera ,fill=(0, 0, 0),font=fonte)
            desenho.text((815,315),z_med ,fill=(0, 0, 0),font=fonte)
            desenho.text((1120,315),exposicao ,fill=(0, 0, 0),font=fonte)
            desenho.text((815,425),lat_med ,fill=(0, 0, 0),font=fonte)
            desenho.text((1120,425),long_med ,fill=(0, 0, 0),font=fonte)
    
            buffer_info = BytesIO()
            infos.save(buffer_info, format="PNG")
            buffer_info.seek(0)
    
            imagem_infos = Image.open(buffer_info)
    
            w_i = w*0.97
            h_i = w_i*alt_info/larg_info
    
            #CRIAR PDF
            pdf = canvas.Canvas(self.pasta+"/fotoindice_"+self.projeto+".pdf", pagesize=(w+64,h+h_i+64))
            x = 32 # Posição X do gráfico no documento PDF)
            y = 32+h_i # Posição Y do gráfico no documento PDF
    
            #ADICIONAR FOTOINDICES
            pdf.drawImage(ImageReader(imagem_pil), x, y, width=w, height=h)
    
            #ADICIONANDO INFORMAÇÕES NO RODAPÉ
            pdf.drawImage(ImageReader(imagem_infos), x+int(w*0.03), y-h_i-2, width=w_i, height=h_i)
    
            pdf.save()
            infos.close()
            imagem_infos.close()
            buffer.close()
            buffer_info.close()
            plt.close()
            imagem_pil.close()

            self.progress_changed.emit(100, '',False, 5)
        except utm.error.OutOfRangeError:
            self.progress_changed.emit(157,'',True,5)
    
        self.etapa_changed.emit(3)

class Thread_6(QThread):
         # Sinais para comunicação entre as threads
         finished = pyqtSignal()
         result_ready = pyqtSignal()
         progress_changed = pyqtSignal(int, str, bool, int)
         etapa_changed = pyqtSignal(int)
     
         def __init__(self, geotiffs, pasta, projeto):
             super().__init__()
             self.geotiffs = geotiffs
             self.pasta = pasta
             self.projeto = projeto

         def run(self):
            # Lógica da thread
            result = self.salvar_geotifs()

            # Emitir o sinal com o resultado
            self.result_ready.emit()

            # Emitir o sinal de finalização
            self.finished.emit()

         def salvar_geotifs(self):
            geotifs = self.geotiffs
            pasta = self.pasta

            if not os.path.exists(pasta+'/Geotiffs'):
                # Cria a pasta se ela não existe
                os.makedirs(pasta+'/Geotiffs')
                pasta = pasta+'/Geotiffs'
            else:
                pasta = pasta+'/Geotiffs'
            for i, geotiff in enumerate(geotifs):
                nome_arquivo_saida = f"{self.projeto}_geotiff_{i}.tif"
                caminho_arquivo_saida = os.path.join(pasta, nome_arquivo_saida)
                
                with open(caminho_arquivo_saida, "wb") as arquivo_saida:
                    arquivo_saida.write(geotiff.read())  # Salva os dados do geotiff no arquivo de saída
                        
                        

  
# Form implementation generated from reading ui file 'qt_fotoindice.ui'
#
# Created by: PyQt6 UI code generator 6.4.2
#
# WARNING: Any manual changes made to this file will be lost when pyuic6 is
# run again.  Do not edit this file unless you know what you are doing.

#diretorio_atual = 'D:/Gerador de Fotoindice/FotoIndice/AeroSAT - Gerador de FotoIndice'
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QApplication, QWidget, QPushButton, QVBoxLayout, QFileDialog
from PyQt6.QtGui import QColor, QPalette, QRegularExpressionValidator
from datetime import datetime
class Ui_Dialog(object):
    def setupUi(self, Dialog):
        regex = QtCore.QRegularExpression("^-?\d+(\.\d+)?$")  # Permite números inteiros ou decimais (exemplo: 10, -5, 3.14, -2.5)
        validator = QtGui.QRegularExpressionValidator(regex)
        Dialog.setObjectName("Dialog")
        Dialog.resize(708, 572)
        self.progresso_funcao = QtWidgets.QProgressBar(parent=Dialog)
        self.progresso_funcao.setGeometry(QtCore.QRect(30, 70, 351, 23))
        self.progresso_funcao.setProperty("value", 0)
        #self.progresso_funcao.setStyleSheet("QProgressBar::chunk { background-color: #00BFFF; width: 1px}" "QProgressBar { text-align: center;}")
        self.progresso_funcao.setObjectName("progresso_funcao")
        self.progresso_total = QtWidgets.QProgressBar(parent=Dialog)
        self.progresso_total.setGeometry(QtCore.QRect(30, 150, 351, 23))
        self.progresso_total.setProperty("format", "%v/%m")
        self.progresso_total.setValue(0)
        self.progresso_total.setRange(0, 3)
        #self.progresso_total.setStyleSheet("QProgressBar::chunk { background-color: #00BFFF;}" "QProgressBar { text-align: center;}")
        self.progresso_total.setObjectName("progresso_total")
        self.label_total = QtWidgets.QLabel(parent=Dialog)
        self.label_total.setGeometry(QtCore.QRect(130,150,150,20))
        self.label_total.setObjectName("label_total")
        self.checkbox = QtWidgets.QCheckBox(parent=Dialog)
        self.checkbox.setChecked(False)  # Define o estado inicial da caixa de checagem
        self.checkbox.setGeometry(QtCore.QRect(30, 180, 15, 15))
        self.checkbox.setObjectName("check_geotif")
        self.check_label = QtWidgets.QLabel(parent=Dialog)
        self.check_label.setGeometry(QtCore.QRect(50,178,80,20))
        self.check_label.setObjectName("check_label")
        self.projeto_edit = QtWidgets.QLineEdit(parent=Dialog)
        self.projeto_edit.setGeometry(QtCore.QRect(510, 30, 131, 20))
        self.projeto_edit.setMaxLength(20)
        self.projeto_edit.setObjectName("projeto_edit")
        self.projeto_label = QtWidgets.QLabel(parent=Dialog)
        self.projeto_label.setGeometry(QtCore.QRect(554, 10, 41, 20))
        self.projeto_label.setObjectName("projeto_label")
        self.data_agora = datetime.now()
        self.data_edit = QtWidgets.QDateEdit(parent=Dialog)
        self.data_edit.setGeometry(QtCore.QRect(440, 70, 121, 21))
        self.data_edit.setDateTime(QtCore.QDateTime(QtCore.QDate(self.data_agora.year, self.data_agora.month, self.data_agora.day), QtCore.QTime(0, 0, 0)))
        self.data_edit.setMinimumDateTime(QtCore.QDateTime(QtCore.QDate(1990, 1, 1), QtCore.QTime(0, 0, 0)))
        self.data_edit.setMinimumDate(QtCore.QDate(1990, 1, 1))
        self.data_edit.setCalendarPopup(True)
        self.data_edit.setObjectName("data_edit")
        self.label = QtWidgets.QLabel(parent=Dialog)
        self.label.setGeometry(QtCore.QRect(470, 50, 61, 16))
        self.label.setObjectName("label")
        self.gsd_edit = QtWidgets.QLineEdit(parent=Dialog)
        self.gsd_edit.setValidator(validator)
        self.gsd_edit.setGeometry(QtCore.QRect(580, 70, 121, 20))
        self.gsd_edit.setObjectName("gsd_edit")
        self.gsd_fuso = QtWidgets.QLabel(parent=Dialog)
        self.gsd_fuso.setGeometry(QtCore.QRect(620, 50, 47, 16))
        self.gsd_fuso.setObjectName("gsd_fuso")
        self.fuso_edit = QtWidgets.QComboBox(parent=Dialog)
        self.fuso_edit.setGeometry(QtCore.QRect(440, 110, 121, 22))
        self.fuso_edit.setLayoutDirection(QtCore.Qt.LayoutDirection.LeftToRight)
        self.fuso_edit.setObjectName("fuso_edit")
        self.fuso_edit.addItem("")
        self.fuso_edit.addItem("")
        self.fuso_edit.addItem("")
        self.fuso_edit.addItem("")
        self.fuso_edit.addItem("")
        self.fuso_edit.addItem("")
        self.fuso_edit.addItem("")
        self.fuso_edit.addItem("")
        self.fuso_labek = QtWidgets.QLabel(parent=Dialog)
        self.fuso_labek.setGeometry(QtCore.QRect(490, 90, 31, 20))
        self.fuso_labek.setObjectName("fuso_labek")
        self.tifs_button = QtWidgets.QPushButton(parent=Dialog)
        self.tifs_button.setGeometry(QtCore.QRect(580, 110, 121, 23))
        self.tifs_button.setAcceptDrops(True)
        self.tifs_button.setObjectName("tifs_button")
        self.csv_button = QtWidgets.QPushButton(parent=Dialog)
        self.csv_button.setGeometry(QtCore.QRect(440, 150, 121, 23))
        self.csv_button.setAcceptDrops(True)
        self.csv_button.setObjectName("csv_button")
        self.kml_button = QtWidgets.QPushButton(parent=Dialog)
        self.kml_button.setGeometry(QtCore.QRect(580, 150, 121, 23))
        self.kml_button.setAcceptDrops(True)
        self.kml_button.setObjectName("kml_button")
        self.gerar_button = QtWidgets.QPushButton(parent=Dialog)
        self.gerar_button.setGeometry(QtCore.QRect(530, 180, 75, 23))
        self.gerar_button.setObjectName("gerar_button")
        self.gerar_button.setEnabled(False)
        self.textBrowser = QtWidgets.QTextBrowser(parent=Dialog)
        self.textBrowser.setGeometry(QtCore.QRect(10, 210, 691, 351))
        self.textBrowser.setObjectName("textBrowser")
        self.kml = ''
        self.counter = 0
        self.pasta = ''
        self.progress_1 = 0
        self.progress_2 = 0
        self.signal_2_3 = pyqtSignal()
        self.counter_thread = 0
        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)
        #habilitar o botão gerar
        # Conectar o sinal textChanged de cada entrada ao método onTextChanged
        self.projeto_edit.textChanged.connect(self.onTextChanged)
        self.gsd_edit.textChanged.connect(self.onTextChanged)
        self.tifs_button.clicked.connect(self.onSelectTIFS)
        self.csv_button.clicked.connect(self.onSelectCSV)
        self.kml_button.clicked.connect(self.onSelectKML)
        self.gerar_button.clicked.connect(self.gerar_imagens_coords)
    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "AeroSAT - Gerador de Fotoindice"))
        self.projeto_label.setText(_translate("Dialog", "Projeto"))
        self.label.setText(_translate("Dialog", "Data de voo"))
        self.gsd_fuso.setText(_translate("Dialog", "GSD (cm)"))
        self.fuso_edit.setItemText(0, _translate("Dialog", "25"))
        self.fuso_edit.setItemText(1, _translate("Dialog", "24"))
        self.fuso_edit.setItemText(2, _translate("Dialog", "23"))
        self.fuso_edit.setItemText(3, _translate("Dialog", "22"))
        self.fuso_edit.setItemText(4, _translate("Dialog", "21"))
        self.fuso_edit.setItemText(5, _translate("Dialog", "20"))
        self.fuso_edit.setItemText(6, _translate("Dialog", "19"))
        self.fuso_edit.setItemText(7, _translate("Dialog", "18"))
        self.fuso_labek.setText(_translate("Dialog", "Fuso"))
        self.tifs_button.setText(_translate("Dialog", "Imagens"))
        self.csv_button.setText(_translate("Dialog", "CP"))
        self.kml_button.setText(_translate("Dialog", "KML"))
        self.gerar_button.setText(_translate("Dialog", "Gerar"))
        self.check_label.setText(_translate("Dialog", "Salvar GEOTIFF"))
    def onTextChanged(self):
        # Verificar se todas as 4 entradas estão preenchidas e o arquivo foi selecionado
        if (
            self.projeto_edit.text() and
            self.gsd_edit.text() and
            self.counter == 2
        ):
            # Habilitar o botão "Enviar" se todas as entradas estiverem preenchidas e o arquivo foi selecionado
            self.gerar_button.setEnabled(True)
        else:
            # Desabilitar o botão "Enviar" se alguma das entradas estiver vazia ou o arquivo não foi selecionado
            self.gerar_button.setEnabled(False)
    def onSelectTIFS(self):
        if not self.pasta:
         # Abrir a caixa de diálogo de seleção de arquivos
             tifs, _ = QFileDialog.getOpenFileNames(Dialog, "Selecionar Imagens", "", "*.tif;*.IIQ")
        else:
             tifs, _ = QFileDialog.getOpenFileNames(Dialog, "Selecionar Imagens", self.pasta, "*.tif;*.IIQ")
         # Exibir o nome do arquivo selecionado na primeira entrada de texto
        if tifs:
            if self.counter == 0 or self.counter == 1:
                self.counter = self.counter + 1
            if not hasattr(Ui_Dialog,'tifs'):
                self.tifs = tifs
                self.pasta, _= os.path.split(tifs[0])
            else:
                self.tifs += tifs
        elif not hasattr(Ui_Dialog,'tifs'):
             if self.counter != 0:
                 self.counter = self.counter - 1
         # Atualizar o status do botão "Enviar"
        self.onTextChanged()
    def onSelectCSV(self):
        if not self.pasta:
            self.csv, _ = QtWidgets.QFileDialog.getOpenFileName(Dialog, 'Selecionar Arquivo dos CP',"*.csv;*.txt")
        else:
            self.csv, _ = QtWidgets.QFileDialog.getOpenFileName(Dialog, 'Selecionar Arquivo dos CP', self.pasta,"*.csv;*.txt")
        if self.csv:
            self.pasta, _= os.path.split(self.csv)
            if self.counter == 0 or self.counter == 1:
               self.counter = self.counter + 1
        else:
            if self.counter != 0:
                self.counter = self.counter - 1
        self.onTextChanged()
    def onSelectKML(self):
        if not self.pasta:
            self.kml, _ = QtWidgets.QFileDialog.getOpenFileName(Dialog, 'Selecionar KML',"*.kml")
        else:
            self.kml, _ = QtWidgets.QFileDialog.getOpenFileName(Dialog, 'Selecionar KML',self.pasta,"*.kml")
    def on_finished_trheads_1(self, result):
        self.imagens_coords, self.size, self.bds, self.escala, self.ind = result
        self.gerar_geotiffs()
    def on_finished_threads_2_3_4(self, result, value):
        if value == 1:
            self.geotiffs_1 = result
        if value == 2:
            self.geotiffs_2 = result
        if value == 3:
            self.geotiffs_3 = result
        self.counter_thread += 1
        if self.counter_thread == 3:
            self.progresso_total.setValue(2)
            self.progresso_funcao.setValue(0)
            self.counter_thread = 0
            self.gerar_fotoindice()
    def gerar_imagens_coords(self):
        self.progresso_total.setValue(0)
        self.projeto = self.projeto_edit.text()
        self.gsd = self.gsd_edit.text()
        self.fuso = self.fuso_edit.currentText()
        self.data = (self.data_edit.date()).toString("dd/MM/yyyy")
        self.progresso_funcao.setValue(0) # PARA CASOS EM QUE SE GERA MAIS DE UMA VEZ SEM EXECUTAR NOVAMENTE O PROGRAMA
        #Etapa 1 - Imagens e cps, tamanho da thumb e extensão do mapa
        self.thread_1 = Thread_1(self.tifs, self.csv, self.gsd, self.fuso)
        self.thread_1.progress_changed.connect(self.on_thread_progress_changed)
        self.thread_1.etapa_changed.connect(self.on_thread_etapa_changed)
        self.thread_1.result_ready.connect(self.on_finished_trheads_1)
        _translate = QtCore.QCoreApplication.translate
        self.label_total.setText(_translate("Dialog", "Verificando CP's das imagens..."))
        self.thread_1.start()
    def gerar_geotiffs(self):
        tif_terco = len(self.tifs) // 3
        tifs_1 = self.tifs[:tif_terco]
        tifs_2 = self.tifs[tif_terco:2*tif_terco]
        tifs_3 = self.tifs[2*tif_terco:]
        self.thread_2 = Thread_2(tifs_1, self.size, self.fuso, int(self.gsd), self.imagens_coords, 1)
        self.thread_2.progress_changed.connect(self.on_thread_progress_changed)
        self.thread_2.etapa_changed.connect(self.on_thread_etapa_changed)
        self.thread_3 = Thread_3(tifs_2, self.size, self.fuso, int(self.gsd), self.imagens_coords, 2)
        self.thread_3.progress_changed.connect(self.on_thread_progress_changed)
        self.thread_3.etapa_changed.connect(self.on_thread_etapa_changed)
        self.thread_4 = Thread_4(tifs_3, self.size, self.fuso, int(self.gsd), self.imagens_coords, 3)
        self.thread_4.progress_changed.connect(self.on_thread_progress_changed)
        self.thread_4.etapa_changed.connect(self.on_thread_etapa_changed)
        self.thread_2.result_ready.connect(self.on_finished_threads_2_3_4)
        self.thread_3.result_ready.connect(self.on_finished_threads_2_3_4)
        self.thread_4.result_ready.connect(self.on_finished_threads_2_3_4)
        self.thread_2.start()
        self.thread_3.start()
        self.thread_4.start()
        _translate = QtCore.QCoreApplication.translate
        self.label_total.setText(_translate("Dialog", "Gerando fotoindices..."))
        #self.signal_2_3.finished.connect(self.gerar_fotoindice)
    def gerar_fotoindice(self):
        self.geotiffs = self.geotiffs_1 + self.geotiffs_2 + self.geotiffs_3
        self.geotiffs.reverse()
        if self.checkbox.isChecked():
           self.salvar_geotiffs()
        self.thread_5 = Thread_5(self.geotiffs, self.kml, self.pasta, self.projeto, self.data, self.gsd, self.fuso, self.escala, 'EPSG:'+str(31960+int(self.fuso)), self.bds, self.ind, diretorio_atual, self.tifs, self.imagens_coords)
        self.thread_5.progress_changed.connect(self.on_thread_progress_changed)
        self.thread_5.etapa_changed.connect(self.on_thread_etapa_changed)
        self.thread_5.start()
        _translate = QtCore.QCoreApplication.translate
        self.label_total.setText(_translate("Dialog", "Criando PDF do fotoindice..."))
        self.tifs=[]
    def salvar_geotiffs(self):
        self.thread_6 = Thread_6(self.geotiffs, self.pasta, self.projeto)
        self.thread_6.progress_changed.connect(self.on_thread_progress_changed)
        self.thread_6.etapa_changed.connect(self.on_thread_etapa_changed)
        self.thread_6.start()
        
        '''
        for i, geot in enumerate(geotiffs):
            geot = geot.open()
            transform = geot.transform
            file_name = self.tifs[i]
            tfw_file_name = file_name.replace('.tif', '.tfw')
            
            with open(tfw_file_name, 'w') as tfw_file:
                tfw_file.write(f"{transform.a:.6f}\n{transform.b:.6f}\n{transform.d:.6f}\n{transform.e:.6f}\n{transform.xoff:.10f}\n{transform.yoff:.10f}")
        '''

    #@pyqtSlot(int)
    def on_thread_progress_changed(self, value, value_str, erro, id):
        cursor = self.textBrowser.textCursor()
        cursor.movePosition(QtGui.QTextCursor.MoveOperation.End)
        # Slot chamado quando o progresso da thread é atualizado
        if id == 2:
            self.progress_1 = value
        elif id == 3:
            self.progress_2 = value
        elif id == 4:
            self.progress_3 = value
        else:
            self.progresso_funcao.setValue(value)
        if id == 5 and value == 100:
                self.textBrowser.setTextColor(QColor("blue"))
                self.textBrowser.insertPlainText('Fotoindice salvo no diretório ' + self.pasta)
                self.textBrowser.setTextCursor(cursor)
                _translate = QtCore.QCoreApplication.translate
                self.label_total.setText(_translate("Dialog", 'Fotoindice salvo no diretório'))
        try:
            if self.thread_2.isRunning() and self.thread_3.isRunning() and self.thread_4.isRunning():
                pgr = (self.progress_1+self.progress_2+self.progress_3)/3
                self.progresso_funcao.setValue(int(pgr))
                if erro:
                    self.textBrowser.setTextColor(QColor("red"))
                    self.textBrowser.insertPlainText('Não foi possível criar o fotoindice da imagem ' + value_str +'\n')
                    self.textBrowser.setTextCursor(cursor)
                else:
                    self.textBrowser.setTextColor(QColor("green"))
                    self.textBrowser.insertPlainText('Fotoindice da imagem ' + value_str + ' gerada com sucesso'+'\n')
                    self.textBrowser.setTextCursor(cursor)

        except:
            pass
        if value == 150:
            self.textBrowser.setTextColor(QColor("red"))
            self.textBrowser.insertPlainText('Há um problema com a formatação do arquivo! Certifique-se que o espaço entre cada valor da tabela seja de apenas um espaço'+'\n')
            self.textBrowser.setTextCursor(cursor)
        if value == 151:
            self.textBrowser.setTextColor(QColor("red"))
            self.textBrowser.insertPlainText('Há um problema com a formatação do arquivo! Verifique se o separador entre cada valor seja um espaço'+'\n')
            self.textBrowser.setTextCursor(cursor)
        if value == 152:
            self.textBrowser.setTextColor(QColor("red"))
            self.textBrowser.insertPlainText('Não foram encontrados CPs para as imagens'+'\n')
            self.textBrowser.setTextCursor(cursor)
        if value==157:
            self.textBrowser.setTextColor(QColor("red"))
            self.textBrowser.insertPlainText('Houve um problema com a extensão da área. Verifique se o fuso indicado é correto.'+'\n')
            self.textBrowser.setTextCursor(cursor)

       
    #@pyqtSlot(int)
    def on_thread_etapa_changed(self, value):
        #self.progresso_total.setValue(value)
        self.progresso_total.setValue(value)

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    app_icon = QtGui.QIcon()
    app_icon.addFile(diretorio_atual + '/resources/aerologo_fotoindice.ico', QtCore.QSize(48,48))
    app.setWindowIcon(app_icon)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec())
