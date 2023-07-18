# Form implementation generated from reading ui file 'gerador_ui_v2.ui'
#
# Created by: PyQt6 UI code generator 6.4.2
#
# WARNING: Any manual changes made to this file will be lost when pyuic6 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QDialog
from PyQt6.QtCore import QThread, pyqtSignal, pyqtSlot, QRegularExpression
from PyQt6.QtGui import QColor, QPalette, QRegularExpressionValidator
import inspect
import os
import geopandas as gpd
import pandas as pd
import shapely as shp
import ee
import utm
import math as mt
import fiona
import json
from fiona.drvsupport import supported_drivers
import gpxpy
import gpxpy.gpx
import simplekml
from shapely.geometry import Point, Polygon
from shapely.affinity import affine_transform
import numpy as np
import threading as thr

#diretorio_atual = os.path.dirname(inspect.getfile(inspect.currentframe()))
#diretorio_atual = os.path.abspath(os.path.dirname(__file__))
#diretorio_atual = os.path.dirname(os.path.abspath(__file__))
diretorio_atual = os.path.abspath(os.path.dirname(sys.executable))

ee_account='' #CONTA OCULTADA
credentials = ee.ServiceAccountCredentials(ee_account, diretorio_atual + '/resources/private-key.json')
ee.Initialize(credentials)



class MyThread(QThread):
    # Sinais para comunicação entre as threads
    finished = pyqtSignal()
    result_ready = pyqtSignal(int)
    progress_changed = pyqtSignal(int)

    def __init__(self, loc_fn, loc, camm, proj, fuso, vel, gsd, soblong, soblat, freq, pot):
        super().__init__()
        self.loc_fn = loc_fn
        self.loc = loc
        self.camm = camm
        self.proj = proj
        self.fuso = fuso
        self.vel = vel
        self.gsd = gsd
        self.soblong = soblong
        self.soblat = soblat
        self.freq = freq
        self.pot = pot

    def run(self):
        # Lógica da thread
        result = self.gerar_pl()

        # Emitir o sinal com o resultado
        self.result_ready.emit(result)

        # Emitir o sinal de finalização
        self.finished.emit()


    def gerar_pl(self):
            def cameraFun(cam):
                class Camera:
                      def __init__(self, f, tp, ccdx, ccdy):
                         self.focal = f #focal em mm
                         self.tp = tp #Tamanho do pixel em mm
                         self.CCDx = ccdx #CCD em X (Pixel)
                         self.CCDy = ccdy  #CCD em Y


                global camera
                if cam == 'PhaseOne 100MP':
                    camera = Camera() #PARAMETROS DAS CAMERAS OCULTADOS
                if cam == 'PhaseOne 80MP':
                    camera = Camera() #PARAMETROS DAS CAMERAS OCULTADOS
                if cam=='Hasselblad':
                    camera = Camera() #PARAMETROS DAS CAMERAS OCULTADOS

                return camera


            class Aeronave:
                 def __init__(self, vel):
                     self.vel = vel #Velocidade da Aeronave em nós



            class Plano:
                 def __init__(self, gsd, longitudinal, lateral, freq, pot):
                     self.gsd = gsd*10 #cm p mm
                     self.long = longitudinal # em %
                     self.lat = lateral # em %
                     self.freq = freq
                     self.pot = pot

            #Detalhe: como os voos da Aerosat são feitos com a câmera em modo paisagem, o recobrimento lateral será maior que o longitudinal porque a cânmera está deitada
            class Resultados:
                def __init__(self, aeronave, camera, plano):
                    self.aero = aeronave
                    self.cam = camera
                    self.plano = plano
                def abrangenciax(self): #Abrangência em X de cada foto
                         return self.cam.CCDx*self.plano.gsd

                def abrangenciay(self): #Abrangência em Y de cada foto
                         return self.cam.CCDy*self.plano.gsd

                def alturaemmetros(self): #Altura do voo (H) em metros
                         return ((self.cam.focal*self.plano.gsd*self.cam.CCDx)/self.cam.tp)/10000000

                def alturaempes(self): #H em pés
                         return self.alturaemmetros()/0.3048

                def recobrimentolat(self): #Recobrimento lateral em metros
                          return ((self.cam.CCDx*self.cam.tp)/self.cam.focal)*self.alturaemmetros()

                def recobrimentolong(self): #Recobrimento longitudinal em metros
                          return ((self.cam.CCDy*self.cam.tp)/self.cam.focal)*self.alturaemmetros()

                def intervalo(self): #Intervalo de tomada das fotos em segundos
                          return (((100-self.plano.long)/100)*self.recobrimentolong())/(self.aero.vel*(1.852/3.6))
                def aerobase(self): #Distância entre as tomadas das fotos em metros
                          return self.recobrimentolong()-self.recobrimentolong()*self.plano.long/100




            def geracp(fxs, base):


                #self.Infos.setText('Aguarde, criando os CPs\n')

                dist = lambda x1,x2,y1,y2: ((x2-x1)**2 + (y2-y1)**2)**(1/2)


                def az(x1,x2,y1,y2):
                    if y2==y1:
                         if x2>x1:
                             az = 90
                         else:
                             az = 270
                    elif x2==x1:
                         if y2>y1:
                             az = 0
                         else:
                             az = 180
                    else:
                        ang = mt.atan((x2-x1)/(y2-y1))*180/mt.pi
                        if x1 < x2 and y1<y2:
                           az = ang
                        elif x1 < x2 and y1>y2:
                           az = 180 + ang
                        elif x1 > x2 and y1>y2:
                           az = 180 + ang
                        else:
                           az = 360 + ang

                    return az

                faixas_cps = []
                for fx in fxs:
                    fxcps = {}
                    cps = []

                    fxcps['Faixa'] = fx['properties']['Layer']

                    x1 = fx['geometry']['coordinates'][0][0]
                    x2 = fx['geometry']['coordinates'][1][0]
                    y1 = fx['geometry']['coordinates'][0][1]
                    y2 = fx['geometry']['coordinates'][1][1]


                    azmt = az(x1,x2,y1,y2)
                    distan = dist(x1,x2,y1,y2)
                    ncp = (round(distan/base))-1 #pq a coordenada inicial e final já são 2 CP's. Numero de CP por faixa

                    fxcps['Azimute'] = azmt
                    fxcps['Número de Fotos'] = ncp
                    fxcps['Distância'] = distan/1000

                    cps.append((x1,y1))
                    for i in range(0,ncp):
                        fator = i+1
                        distancia = base*fator

                        x = x1 + distancia*mt.sin(azmt*mt.pi/180)
                        y = y1 + distancia*mt.cos(azmt*mt.pi/180)

                        cps.append((x,y))


                    cps.append((x2,y2))




                    fxcps['CPs'] = cps
                    faixas_cps.append(fxcps)

                    #self.Infos.setText('Os CPs foram criados\n')



                self.progress_changed.emit(16)
                return faixas_cps

            def getaltitude(faixas, fuso):

                    #self.Infos.setText('Aguarde, estamos coletando as altitudes médias de cada faixa\n')


                    def lineToPoints(linha, cont):
                      length = linha.length()
                      step = length.divide(cont)
                      dist = ee.List.sequence(0, length, step)

                      def tomap(s):
                          line = ee.List(s).get(0)
                          offset = ee.List(s).get(1)
                          return makePointFeature(ee.Geometry(line).coordinates().get(0), offset)

                      def makePointFeature(coord, offset):
                          pt = ee.Algorithms.GeometryConstructors.Point(coord)
                          coords = pt.coordinates()
                          return ee.Feature(pt).set({'offset': offset, 'lat': coords.get(0), 'lon': coords.get(1)})

                      lines = linha.cutLines(dist).geometries()
                      points = lines.zip(dist).map(tomap)
                      points = points.add(makePointFeature(linha.coordinates().get(-1), length))

                      return ee.FeatureCollection(points)


                    def images():
                         #Imagens
                         imagem_2 = ee.Image('CGIAR/SRTM90_V4')
                         elevacao = imagem_2.select('elevation')
                         slope = ee.Terrain.slope(elevacao)
                         topo = ee.Image.cat(elevacao, slope)
                         topocol = ee.ImageCollection([topo])

                         return topocol


                    def geom_csv(faixas, zona):
                            def rasterExtraction(image):
                              try:
                                 feature = image.sampleRegions(collection = points, geometries = False)
                                 return feature
                              except:
                                  self.progress_changed.emit(404)
                                  sys.exit()


                            topocol = images()
                            #gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
                            #df = gpd.read_file('C:/Users/User/Downloads/gee/teste.kml', driver='KML')
                            #df_json = df.to_json()
                            valores = []
                            for faixa in faixas:
                               if len(faixa['geometry']['coordinates']) == 2:
                                     p1=faixa['geometry']['coordinates'][0] #Pensar em transformar as coordenadas antes p/ wgs84
                                     p2=faixa['geometry']['coordinates'][1]
                                     p1_wgs84=list(utm.to_latlon(p1[0],p1[1],zona,northern=False))
                                     p2_wgs84=list(utm.to_latlon(p2[0],p2[1],zona,northern=False))

                                     linha = ee.Geometry.LineString([[p1_wgs84[1],p1_wgs84[0]],[p2_wgs84[1],p2_wgs84[0]]])
                                     points = lineToPoints(linha, int(linha.length().getInfo())*4/1000)

                                     results = topocol.filterBounds(points).map(rasterExtraction).flatten()
                                     resultados = results.getInfo()


                                     valores.append({'Faixa': faixa['properties']['Layer'], 'Features': resultados['features']})




                            data = []
                            for valor in valores:
                                  elev = []
                                  feats = valor['Features']
                                  for feat in feats:

                                          #data.append([valor['Faixa'],feat['properties']['elevation']])
                                      elev.append(feat['properties']['elevation'])

                                  try:
                                      data.append([valor['Faixa'],sum(elev)/len(elev)])
                                  except ZeroDivisionError:
                                      data.append([valor['Faixa'],0])

                                  #self.Infos.setText('A altitude média da da '+valor['Faixa']+' é '+str(sum(elev)/len(elev)))



                            return data

                    #self.Infos.setText('As altitudes foram coletadas\n')


                    self.progress_changed.emit(32)
                    return geom_csv(faixas, fuso)


            def create_rout(linestrings,zona,output, proj):
                #self.Infos.setText('Aguarde, estamos criando as rotas')
                # converter as linestrings para uma rota GPX
                gpx = gpxpy.gpx.GPX()
                for i, row in linestrings.iterrows():
                    route = gpxpy.gpx.GPXRoute()
                    route.name = linestrings['Layer'][i]
                    for point in row['geometry'].coords:
                        coord_wgs84=list(utm.to_latlon(point[0],point[1],zona,northern=False))
                        route.points.append(gpxpy.gpx.GPXRoutePoint(latitude=coord_wgs84[0], longitude=coord_wgs84[1]))
                    gpx.routes.append(route)

                # salvar o arquivo GPX
                with open(output+'/rota_'+proj+'.gpx', 'w') as f:
                    f.write(gpx.to_xml())

                #self.Infos.setText('As rotas foram criadas\n')

                self.progress_changed.emit(64)




            def generate_rectangles(cps, length, width, output, epsg,proj):
                #self.Infos.setText('Aguarde, estamos criando o recobrimento\n')
                buffer_list = []
                cp_list = []

                for cp in cps:

                    point_list = cp['CPs']
                    azimuth = cp['Azimute']

                    azimuth_rad = mt.radians(azimuth)

                    for point in point_list:

                           center = Point(point)
                           rotation_matrix = np.array([[np.cos(azimuth_rad), -np.sin(azimuth_rad)], [np.sin(azimuth_rad), np.cos(azimuth_rad)]])
                           translate_matrix = np.array([center.x, center.y])
                           vertices = np.array([[-length/2, -width/2], [length/2, -width/2], [length/2, width/2], [-length/2, width/2]])
                           vertices = vertices.dot(rotation_matrix) + translate_matrix
                           rectangle = Polygon(vertices)
                           buffer_list.append(rectangle)
                           buffer_list.append(center)


                buffer_gdf = gpd.GeoDataFrame(geometry=buffer_list,crs=epsg)
                fiona.supported_drivers['KML'] = 'rw'
                buffer_gdf.to_file(output+'/recobrimento_'+proj+'.kml',driver='KML')

                #self.Infos.setText('O recobrimento foi criado\n')
                self.progress_changed.emit(75)



            def cria_csv(res, cp, altitudes, proj, loc, camm,fuso):
                 #self.Infos.setText('Aguarde, estamos criando o relatório\n')
                 dados_gerais = [['','Geral',''],
                                ['Projeto',proj],
                                ['Camera',camm],
                                ['GSD (cm)', "{:.1f}".format(res.plano.gsd/10).replace('.', ',')],
                                ['Sobreposicao Longitudinal (%)', int(res.plano.long)],
                                ['Sobreposicao Lateral (%)',int(res.plano.lat)],
                                ['Intervalo de Tomada (s)', "{:.3f}".format(res.intervalo()).replace('.', ',')],
                                [''],
                                ['Frequencia Laser (kHz)', int(plano.freq)],
                                ['Potencia Laser (%)', int(plano.pot)],
                                [''],
                                ['Altura de Voo (pes)', int(round(res.alturaemmetros()*3.28))],
                                ['Altura de Voo (m)',"{:.3f}".format(res.alturaemmetros()).replace('.', ',')],
                                ['Velocidade (nos)', int(aeronave.vel)],
                                [''],
                                [''],
                                ['Faixa', 'Altitude de Voo (pes)']]

                 altura = res.alturaemmetros()
                 temp = []
                 for i in range(0,len(altitudes)):
                        temp.append([altitudes[i][0], int((altitudes[i][1]+altura)*3.28)])


                 temp.sort()
                 for fxx in temp:
                     dados_gerais.append(fxx)


                 csv = pd.DataFrame(dados_gerais)
                 csv.to_csv(loc+'/relatorio_'+proj+'.csv',sep=';',index=False,header=False,encoding='UTF-8',decimal=',',float_format='%.3f')
                 #self.Infos.setText('O relatório foi gerado\n')
                 self.progress_changed.emit(90)


            def gera_kml_riacquire(gdf,zona,proj,loc):
                        #self.Infos.setText('Aguarde, estamos criando o KMl para o RiACQUIRE\n')

                        # Criar um objeto KML
                        kml = simplekml.Kml()



                        # Para cada linha no GeoDataFrame
                        for idx, row in gdf.iterrows():
                            # Obter o nome do layer da linha
                            layer_name = row["Layer"]
                            # Obter as coordenadas da linha
                            coords = list(row["geometry"].coords)

                            p1 = coords[0]
                            p2 = coords[1]


                            p1_wgs84=list(utm.to_latlon(p1[0],p1[1],zona,northern=False))
                            p2_wgs84=list(utm.to_latlon(p2[0],p2[1],zona,northern=False))

                            coords = [[p1_wgs84[1],p1_wgs84[0]],[p2_wgs84[1],p2_wgs84[0]]]
                            coords_point = [((p1_wgs84[1]+p2_wgs84[1])/2,(p1_wgs84[0]+p2_wgs84[0])/2)]




                            # Criar um objeto de linha no KML com as coordenadas e o nome do layer
                            line = kml.newlinestring(name=layer_name, coords=coords)
                            point = kml.newpoint(name=layer_name, coords=coords_point)

                        # Salvar o arquivo KML
                        kml.save(loc+'/FAIXAS_RCQ_'+proj+".kml")

                        #self.Infos.setText('Todo o plano foi gerado e está salvo em ' + loc+'\n')
                        self.progress_changed.emit(100)



            aeronave = Aeronave(self.vel)
            plano = Plano(self.gsd, self.soblong, self.soblat, self.freq, self.pot)
            cameraa = cameraFun(self.camm)
            result = Resultados(aeronave, cameraa, plano)

            EPSG = lambda fus: 'EPSG:' + str(31978+fus-18)


            faixas_dxf = gpd.read_file(self.loc_fn)
            faixas_gjson = json.loads((faixas_dxf.loc[faixas_dxf.geom_type == 'LineString']).to_json())
            faixas = faixas_gjson['features']
            faixas_dxf = faixas_dxf[faixas_dxf.geometry.type == "LineString"]

            cp = geracp(faixas, result.aerobase())
            altitudes = getaltitude(faixas, self.fuso)
            create_rout(faixas_dxf,self.fuso,self.loc, self.proj)
            generate_rectangles(cp, result.recobrimentolat(), result.recobrimentolong(), self.loc,EPSG(self.fuso), self.proj)
            cria_csv(result, cp, altitudes, self.proj, self.loc, self.camm, self.fuso)
            gera_kml_riacquire(faixas_dxf,self.fuso,self.proj,self.loc)

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(487, 422)
        Dialog.setMinimumSize(QtCore.QSize(487, 422))
        Dialog.setMaximumSize(QtCore.QSize(487, 422))
        Dialog.setStyleSheet("")
        # Definir uma expressão regular para validar a entrada (números inteiros ou decimais)
        regex = QRegularExpression("^-?\d+(\.\d+)?$")  # Permite números inteiros ou decimais (exemplo: 10, -5, 3.14, -2.5)
        validator = QRegularExpressionValidator(regex)
        self.Infos = QtWidgets.QTextBrowser(parent=Dialog)
        self.Infos.setGeometry(QtCore.QRect(50, 380, 391, 31))
        self.Infos.setObjectName("Infos")
        self.Proj = QtWidgets.QLineEdit(parent=Dialog)
        self.Proj.setGeometry(QtCore.QRect(160, 30, 171, 20))
        self.Proj.setStyleSheet("bacground-color: rgb(255,255,255)")
        self.Proj.setObjectName("Proj")
        self.label = QtWidgets.QLabel(parent=Dialog)
        self.label.setGeometry(QtCore.QRect(230, 10, 47, 13))
        self.label.setObjectName("label")
        self.dxf = QtWidgets.QPushButton(parent=Dialog)
        self.dxf.setGeometry(QtCore.QRect(60, 80, 171, 23))
        self.dxf.setAcceptDrops(True)
        self.dxf.setObjectName("dxf")
        self.dxf.clicked.connect(self.botao_dxf_click)
        self.loc = QtWidgets.QPushButton(parent=Dialog)
        self.loc.setGeometry(QtCore.QRect(60, 120, 171, 23))
        self.loc.setObjectName("loc")
        self.loc.clicked.connect(self.botao_loc_click)
        self.camm = QtWidgets.QComboBox(parent=Dialog)
        self.camm.setGeometry(QtCore.QRect(250, 80, 171, 22))
        self.camm.setObjectName("camm")
        self.camm.addItem("")
        self.camm.addItem("")
        self.camm.addItem("")
        self.label_2 = QtWidgets.QLabel(parent=Dialog)
        self.label_2.setGeometry(QtCore.QRect(310, 60, 51, 16))
        self.label_2.setObjectName("label_2")
        self.fuso = QtWidgets.QComboBox(parent=Dialog)
        self.fuso.setGeometry(QtCore.QRect(250, 120, 171, 22))
        self.fuso.setObjectName("fuso")
        self.fuso.addItem("")
        self.fuso.addItem("")
        self.fuso.addItem("")
        self.fuso.addItem("")
        self.fuso.addItem("")
        self.fuso.addItem("")
        self.fuso.addItem("")
        self.fuso.addItem("")
        self.label_3 = QtWidgets.QLabel(parent=Dialog)
        self.label_3.setGeometry(QtCore.QRect(320, 100, 31, 16))
        self.label_3.setObjectName("label_3")
        self.vel = QtWidgets.QLineEdit(parent=Dialog)
        self.vel.setGeometry(QtCore.QRect(310, 190, 131, 20))
        self.vel.setInputMethodHints(QtCore.Qt.InputMethodHint.ImhPreferNumbers)
        self.vel.setObjectName("vel")
        self.vel.setValidator(validator)
        self.gsd = QtWidgets.QLineEdit(parent=Dialog)
        self.gsd.setGeometry(QtCore.QRect(50, 190, 131, 20))
        self.gsd.setInputMethodHints(QtCore.Qt.InputMethodHint.ImhPreferNumbers)
        self.gsd.setObjectName("gsd")
        self.gsd.setValidator(validator)
        self.sublong = QtWidgets.QLineEdit(parent=Dialog)
        self.sublong.setGeometry(QtCore.QRect(50, 240, 131, 20))
        self.sublong.setInputMethodHints(QtCore.Qt.InputMethodHint.ImhPreferNumbers)
        self.sublong.setObjectName("sublong")
        self.sublong.setValidator(validator)
        self.sublat = QtWidgets.QLineEdit(parent=Dialog)
        self.sublat.setGeometry(QtCore.QRect(310, 240, 131, 20))
        self.sublat.setInputMethodHints(QtCore.Qt.InputMethodHint.ImhPreferNumbers)
        self.sublat.setObjectName("sublat")
        self.sublat.setValidator(validator)
        self.label_4 = QtWidgets.QLabel(parent=Dialog)
        self.label_4.setGeometry(QtCore.QRect(40, 170, 141, 20))
        font = QtGui.QFont()
        font.setPointSize(7)
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.label_5 = QtWidgets.QLabel(parent=Dialog)
        self.label_5.setGeometry(QtCore.QRect(300, 170, 141, 20))
        font = QtGui.QFont()
        font.setPointSize(7)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.label_6 = QtWidgets.QLabel(parent=Dialog)
        self.label_6.setGeometry(QtCore.QRect(40, 220, 141, 20))
        font = QtGui.QFont()
        font.setPointSize(7)
        self.label_6.setFont(font)
        self.label_6.setObjectName("label_6")
        self.label_7 = QtWidgets.QLabel(parent=Dialog)
        self.label_7.setGeometry(QtCore.QRect(300, 220, 161, 20))
        font = QtGui.QFont()
        font.setPointSize(7)
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.gerar = QtWidgets.QPushButton(parent=Dialog)
        self.gerar.setGeometry(QtCore.QRect(190, 320, 111, 23))
        self.gerar.setCursor(QtGui.QCursor(QtCore.Qt.CursorShape.PointingHandCursor))
        self.gerar.setStyleSheet("background-color: #3BA0ED;border-width: 2px; border-style: solid; border-color: #000000;border-radius: 5px;color: #000000; font-family: \"Andale Mono\"; font-size: 12px; cursor: pointer;")
        self.gerar.setObjectName("gerar")
        self.gerar.clicked.connect(self.gerar_plano)
        self.freq = QtWidgets.QLineEdit(parent=Dialog)
        self.freq.setGeometry(QtCore.QRect(50, 290, 131, 20))
        self.freq.setInputMethodHints(QtCore.Qt.InputMethodHint.ImhPreferNumbers)
        self.freq.setObjectName("freq")
        self.freq.setValidator(validator)
        self.poten = QtWidgets.QLineEdit(parent=Dialog)
        self.poten.setGeometry(QtCore.QRect(310, 290, 131, 20))
        self.poten.setInputMethodHints(QtCore.Qt.InputMethodHint.ImhPreferNumbers)
        self.poten.setObjectName("poten")
        self.poten.setValidator(validator)
        self.label_8 = QtWidgets.QLabel(parent=Dialog)
        self.label_8.setGeometry(QtCore.QRect(70, 270, 111, 20))
        font = QtGui.QFont()
        font.setPointSize(7)
        self.label_8.setFont(font)
        self.label_8.setObjectName("label_8")
        self.label_9 = QtWidgets.QLabel(parent=Dialog)
        self.label_9.setGeometry(QtCore.QRect(330, 270, 101, 20))
        font = QtGui.QFont()
        font.setPointSize(7)
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")
        self.progressBar = QtWidgets.QProgressBar(parent=Dialog)
        self.progressBar.setGeometry(QtCore.QRect(50, 350, 421, 23))
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName("progressBar")

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "AeroSAT - Planejamento de Voo"))
        self.label.setText(_translate("Dialog", "Projeto"))
        self.dxf.setText(_translate("Dialog", "Selecione o DXF"))
        self.loc.setText(_translate("Dialog", "Salvar..."))
        self.camm.setItemText(0, _translate("Dialog", "PhaseOne 100MP"))
        self.camm.setItemText(1, _translate("Dialog", "PhaseOne 80MP"))
        self.camm.setItemText(2, _translate("Dialog", "Hasselblad"))
        self.label_2.setText(_translate("Dialog", "Câmera"))
        self.fuso.setItemText(0, _translate("Dialog", "18"))
        self.fuso.setItemText(1, _translate("Dialog", "19"))
        self.fuso.setItemText(2, _translate("Dialog", "20"))
        self.fuso.setItemText(3, _translate("Dialog", "21"))
        self.fuso.setItemText(4, _translate("Dialog", "22"))
        self.fuso.setItemText(5, _translate("Dialog", "23"))
        self.fuso.setItemText(6, _translate("Dialog", "24"))
        self.fuso.setItemText(7, _translate("Dialog", "25"))
        self.label_3.setText(_translate("Dialog", "Fuso"))
        self.label_4.setText(_translate("Dialog", "    Velocidade da Aeronave (nós)"))
        self.label_5.setText(_translate("Dialog", "     Ground Sample Distance (cm)"))
        self.label_6.setText(_translate("Dialog", "    Sobreposição Longitudinal (%)"))
        self.label_7.setText(_translate("Dialog", "        Sobreposição Lateral (%)"))
        self.gerar.setText(_translate("Dialog", "Gerar Plano"))
        self.label_8.setText(_translate("Dialog", "Frequência Laser (Hz)"))
        self.label_9.setText(_translate("Dialog", "Potência Laser (%)"))

    def botao_dxf_click(self):
        #global filename_dxf
        #filename_dxf = filedialog.askopenfilename(initialdir = "/",
        #                            title = "Select a File",
        #                            filetypes = (("Arquivo AutoCAD",
        #                                            "*.dxf*"),
        #                                        ("All Files",
        #                                            "*.*")))

        # Change label contents
        self.fileName = QtWidgets.QFileDialog.getOpenFileName(Dialog, 'OpenFile',"*.dxf")
        if not self.fileName:
            self.Infos.setText("Você ainda não selecionou o arquivo DXF\n")
        else:
            self.Infos.setText("Você selecionou o arquivo " + self.fileName[0]+'\n')
        #self.myTextBox.setText(fileName)
        #print(fileName)
    def botao_loc_click(self):
        self.saveName = QtWidgets.QFileDialog.getExistingDirectory(Dialog, 'SaveFile')
        if not self.saveName:
            self.Infos.setText("Você ainda não selecionou a pasta onde os arquivos serão salvos"+'\n')
        else:
            self.Infos.setText("Os arquivos serão salvos em " + self.saveName+'\n')
        #print(saveName)

    def gerar_plano(self):
            loc_fn = self.fileName[0]
            loc = self.saveName
            camm = self.camm.currentText()
            proj = self.Proj.text()
            fuso=int(self.fuso.currentText())
            vel = int(self.gsd.text())
            gsd = float(self.vel.text())
            soblong = int(self.sublong.text())
            soblat = int(self.sublat.text())
            freq = int(self.freq.text())
            pot = int(self.poten.text())
            #(self, loc_fn, loc, camm, proj, fuso, vel, gsd, soblong, soblat, freq, pot):
            self.thread = MyThread(loc_fn, loc, camm, proj, fuso, vel, gsd, soblong, soblat, freq, pot)
            self.thread.progress_changed.connect(self.on_thread_progress_changed)

            self.thread.start()

            self.Infos.setText("Gerando CP's, aguarde...'")


    #@pyqtSlot(int)
    def on_thread_progress_changed(self, value):
        # Slot chamado quando o progresso da thread é atualizado
        if value!=404:
              self.progressBar.setValue(value)
        if value == 16:
            self.Infos.setText("Gerando CP's, aguarde...'")
        if value == 32:
            self.Infos.setText("Gerando Altitudes, aguarde...")
        if value == 64:
            self.Infos.setText("Gerando o arquivo de rota, aguarde...")
        if value == 75:
            self.Infos.setText("Gerando o recobrimento das fotos, aguarde...")
        if value == 90:
            self.Infos.setText("Gerando o relatório, aguarde...")
        #if value == 90:
            #self.Infos.setText("Gerando o kml para o Riacquare, aguarde...")
        if value == 100:
            self.Infos.setText("Todo o plano foi gerado e guardado!")
        if value == 404:
            #palette = QPalette()
            self.progressBar.setValue(0)
            self.Infos.setTextColor(QColor("red"))
            self.Infos.setText("Erro na obtenção das altitudes. Verifique se o fuso setado está correto.")






if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    app_icon = QtGui.QIcon()
    app_icon.addFile(diretorio_atual + '/resources/aerologo_transp_pk6_icon.ico', QtCore.QSize(48,48))
    app.setWindowIcon(app_icon)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec())

    #Botao Select File


