import rasterio as rst
from PIL import Image
import os
import numpy as np
import math as mt
from rasterio.transform import Affine
import csv
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
import matplotlib.ticker as ticker
import time

def coords_cantos(coords, altura, largura, gsd,altura_o,largura_o):
    x = float(coords[0])
    y = float(coords[1])
    ang = float(coords[2])

    dist = mt.sqrt(((altura_o)*gsd)**2 + ((largura_o)*gsd)**2)

    gcp1 = rst.control.GroundControlPoint(0, 0, x + dist*mt.sin((315+ang)*mt.pi/180), y + dist*mt.cos((315+ang)*mt.pi/180))
    gcp2 = rst.control.GroundControlPoint(0, largura, x + dist*mt.sin((ang+45)*mt.pi/180), y + dist*mt.cos((45+ang)*mt.pi/180))
    gcp3 = rst.control.GroundControlPoint(altura/2, largura/2, x, y)
    gcp4 = rst.control.GroundControlPoint(altura, largura, x + dist*mt.sin((ang+135)*mt.pi/180), y + dist*mt.cos((135+ang)*mt.pi/180))
    gcp5 = rst.control.GroundControlPoint(altura, 0, x + dist*mt.sin((ang+225)*mt.pi/180), y + dist*mt.cos((225+ang)*mt.pi/180))

    return [gcp1, gcp2, gcp3, gcp4, gcp5]

def mem_file(tifs,size,imagens_coords,epsg,pixel_size,og_h,og_w):

    geotiffs = []
    for tif in tifs:
        img = Image.open(tif)

        img.thumbnail((size[0],size[1]))

        img_np = np.array(img)

        names = tif.split('\\')

        try:
           xy = imagens_coords[names[1]]
           altura, largura = img_np.shape[:2]
           num_bandas = img_np.shape[2]

           #transform = Affine.translation(float(xy[0]) - (size[0]*pixel_size/2), float(xy[1]) + (size[1]*pixel_size/2))
           #transform.rotation(float(xy[2]))
           gcps = coords_cantos(xy, altura, largura, pixel_size, og_h, og_w)
           transform = rst.transform.from_gcps(gcps)
           img_np = np.transpose(img_np, (2, 0, 1))


           img.close()

           memfile = rst.io.MemoryFile()
           with memfile.open(driver='GTiff',width=largura, height=altura, count=num_bandas, dtype=img_np.dtype,transform=transform,crs=epsg) as dst:
               dst.write(img_np)
           geotiffs.append(memfile)
        except:
            pass
    return geotiffs

def fig_fotoindice(bds,geotiffs,pasta):
	fig, ax = plt.subplots()
    #fig, ax= plt.subplots(figsize=(r1, r2))
    for geotiff_path in geotiffs:
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

        ax.imshow(imagem_rotacionada,extent=[extent.left,extent.right,extent.bottom,extent.top],alpha=imagem_rotacionada[:,:,3]/255)

    formatter = ticker.StrMethodFormatter('{x:.0f}')
    plt.gca().xaxis.set_major_formatter(formatter)
    plt.gca().yaxis.set_major_formatter(formatter)

    plt.xlabel('Coordenada E')
    plt.ylabel('Coordenada N')
    ax.set_xlim(bds[1], bds[0])
    plt.xticks(fontsize=8)
    ax.set_ylim(bds[3], bds[2])
    plt.yticks(fontsize=8,rotation=90,verticalalignment='center')
    plt.grid(linewidth=0.1,color='black')
    #plt.show()
    plt.savefig(pasta+'/fotoindice.png', dpi=1500)




#pasta = 'D:/Gerador de Fotoindice/teste'
pasta = ''
caminhos = [os.path.join(pasta, nome) for nome in os.listdir(pasta)]
arquivos = [arq for arq in caminhos if os.path.isfile(arq)]
tifs = [arq for arq in arquivos if arq.lower().endswith(".tif")]

xn = []
yn = []
with open(pasta+'','r') as file:
    reader = csv.reader(file)
    imagens_coords = {}
    for row in reader:
        try:
            infos = row[0].split(' ')
            if len(infos)>=7:
              xn.append(infos[1])
              yn.append(infos[2])
              imagens_coords[infos[0]]=[infos[1],infos[2], infos[6]]
            else:
              print('Há um problema com a formatação do CSV! Certifique-se que o espaço entre cada valor da tabela seja de apenas um espaço')
              break
        except:
            print('Há um problema com a formatação do CSV! Verifique se o separador entre cada valor seja um espaço')
            break


if len(tifs)!=len(imagens_coords):
    tempo = {}
    xn = []
    yn = []
    for tf in tifs:
        nm = tf.split('\\')[1]
        try:
           vet = imagens_coords[nm]
           tempo[nm]=vet
           xn.append(vet[0])
           yn.append(vet[1])
        except:
               pass
    imagens_coords = {}
    imagens_coords = tempo



if len(imagens_coords)>0:

    epsg = 'EPSG:31983'
    pixel_size = 10/100
    sizes = [float(max(xn)),float(min(xn)),float(max(yn)),float(min(yn))]
    dif_x = int(sizes[0] - sizes[1])*1000
    dif_y = int(sizes[2] - sizes[3])*1000

    img_temp = Image.open(tifs[0])
    og_w = img_temp.width
    og_h = img_temp.height

    if dif_x >= dif_y:
        escala = dif_x/380 # pensando numa folha A3 deitada
        ex_x = int((dif_x/escala)/(pixel_size*10)) #multiplicando por 10 o gsd para que seja o tamanho da foto no A3
        fator = og_w/ex_x
        ex_y = int(og_h/fator)
    else:
        escala = dif_y/350 #pensando numa folha A3 em pé
        ex_y = int((dif_y/escala)/(pixel_size*10))
        fator = og_h/ex_y
        ex_x = int(og_w/fator)

    size = [ex_x,ex_y]


    geotiffs = mem_file(tifs,size,imagens_coords,epsg,pixel_size,og_h,og_w)
    if len(geotiffs)>0:
        bds = [float(max(xn))+2000, float(min(xn))-2000, float(max(yn))+2000, float(min(yn))-2000]
        fig_fotoindice(bds,geotiffs,pasta)
        print(str(time.localtime().tm_hour) + ':' + str(time.localtime().tm_min))
    else:
        print('As imagens não foram georreferenciadas porque não foram encontrados CPs ou algum outro erro ocorreu')
else:
    print('Não foram encontrados CPs para as imagens')
