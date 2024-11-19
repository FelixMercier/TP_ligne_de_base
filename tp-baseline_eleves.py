# -*- coding: utf-8 -*-
"""
TP Compensation de lignes de base GNSS
PPMD 2024

@author: Xavier Collilieux & Julien Barnéoud
"""

import numpy as np
import pandas as pd


class Base:
    """
    Class "Base" contenant un objet ligne de base contenant le vecteur et les covarinces
    """
    
    def __init__(self,code1, code2, X2mX1, cov):
        """
        Constructeur. Liste des attributs
        
        Parameters
        ----------
        code1 : str
            code sta 1
        code2 : str
            code sta 2
        X2mX1 : array (3,)
            X2 - X1
        cov : array (3,3)
            sigx^2, sigxy, sigxz, sigy^2, sigyz, sigz^2 

        Returns
        -------
        None.

        """
        self.code1 = code1
        self.code2 = code2
        self.X2mX1 = X2mX1
        self.cov   = np.zeros((3,3))
        self.cov[0,0] = cov[0]
        self.cov[1,1] = cov[3]
        self.cov[2,2] = cov[5]
        self.cov[1,0] = cov[1]
        self.cov[0,1] = cov[1]
        self.cov[2,0] = cov[2]
        self.cov[0,2] = cov[2]
        self.cov[2,1] = cov[4]
        self.cov[1,2] = cov[4]
    
     
    def affiche(self):
        """ Affichage d'un élément de la classe """
        print('------------------------------------------------------')        
        print('Ligne de base : ')
        print("  Code 1, code2 : %s %s" % (self.code1, self.code2))
        print(self.X2mX1)
        print(self.cov)
        print('------------------------------------------------------')



def lecture_covar(nomfic):
    """
    Lire un fichier ligne de base et le formater sous la forme d'une liste d'objet Base.
    Format fichier :
        # Lecture d'un fichier de ligne de base type LGO
        #
        # MLXC         CTA2            -4206.3904    12289.6029     2539.9962
        # 26.661961940862E-08     17.774641293908E-09     14.219713035126E-08
        # 71.098565175632E-09     17.774641293908E-09
        # 23.107033682080E-08

    Parameters
    ----------
    nomfic : str
        nom du fichier

    Returns
    -------
    listb : list of Base object
        liste contenant des lignes de base (type base) lus dans le fichier

    """
    
    listb = []
    # Ouverture du fichier
    fic  = open(nomfic,'r')
    lines= fic.readlines()
    
    for i in np.arange(0,len(lines),4):
        code1 = lines[i].split()[0]
        code2 = lines[i].split()[1]
        dX    = np.array([float(j) for j in lines[i].split()[2:5]])
        cov   = np.zeros((6,1))        
        cov[0:3]   = np.array([float(j) for j in lines[i+1].split()])[0:3].reshape(3,1)        
        cov[3:5]   = np.array([float(j) for j in lines[i+2].split()])[0:2].reshape(2,1)        
        cov[5]     = np.array(float(lines[i+3].split()[0]))        
        

        bas  = Base(code1, code2, dX, cov)        
        
        # Ajout de l'objet à la liste
        listb.append(bas)
        
    fic.close()
    
    return listb

if __name__ == "__main__":

    # Lecture du fichier
    listB = lecture_covar('lignesdebase.txt')
    print(listB[1].code2)
    # Affichage du premier élément    
    listB[1].affiche()
        
coord_appro = np.genfromtxt('coordappro.txt', delimiter=";", dtype=str)
stat = list(coord_appro[:, 0])
stat.remove('BANON')
stat.remove('VLX1')
stat = np.array(stat)
n = len(stat)

n_lb = len(listB)

B = np.zeros((3*n_lb, 1))
P = np.zeros((3*n_lb, 3*n_lb))
A = np.zeros((3*n_lb, 3*n))
V = np.zeros((3*(n-2), 1))

#Remplissage de A, B, P
for k in range(len(listB)):
    
    #P
    P[3*k:3*k+3, 3*k:3*k+3] = listB[k].cov
    
    #A, B
    stat_connues = ['BANON', 'VLX1']
    #cas1 (connue -> unk)
    if listB[k].code1 in stat_connues :
        j = np.where(stat == listB[k].code2)[0][0]
        A[ 3*k:3*k+3, 3*j:3*j+3 ] = np.eye(3)
        line = np.where(coord_appro == listB[k].code1)[0][0]
        coord = np.zeros((3,1))
        coord[0, 0], coord[1, 0], coord[2, 0] = float(coord_appro[line , 1]), float(coord_appro[line , 2]), float(coord_appro[line , 3])
        B[ 3*k:3*k+3 ] = listB[k].X2mX1.reshape((3,1)) + coord
        
    #cas 2 (unk -> connue) 
    elif listB[k].code2 in stat_connues :
        i = np.where(stat == listB[k].code1)[0][0]
        A[ 3*k:3*k+3, 3*i:3*i+3 ] = -1*np.eye(3)
        line = np.where(coord_appro == listB[k].code2)[0][0]
        coord = np.zeros((3,1))
        coord[0, 0], coord[1, 0], coord[2, 0] = float(coord_appro[line , 1]), float(coord_appro[line , 2]), float(coord_appro[line , 3])
        B[ 3*k:3*k+3 ] = listB[k].X2mX1.reshape((3,1)) - coord
        
        
    #cas 3 (unk -> unk)
    else :
        i, j = np.where(stat == listB[k].code1)[0][0], np.where(stat == listB[k].code2)[0][0]
        A[ 3*k:3*k+3, 3*i:3*i+3 ] = -1*np.eye(3)
        A[ 3*k:3*k+3, 3*j:3*j+3 ] = np.eye(3)


N = A.T @ P @ A
K = A.T @ P @ B
Xchap = np.linalg.inv(N) @ K



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    