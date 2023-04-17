

####################################################################
# @Protein libs
#
#
# -> pmath.py
#        float get_dist(atom1(x,y,z),atom2(x,y,z))                                    # get eucliean distance from 2 atoms
#        list  get_point(vector(x,y,z),ref_point(x,y,z),lambda_v)                     # get a point on a vector from colinear vector, referent point
#                                                                                     - and distance to referent point
#        list   get_vector(pa,pb)                                                      # get colinear vector between 2 points
#        float  get_min_dist(dico1,dico2) dico[""]= [x,y,z]                            # get minimal distance between 2 series of atoms
#        list   get_barycenter(coord) -> [x,y,z]                                       # Give the barycenter of an ensemble of coordinates given in list coord=[[x1,y1,z1],[x2,y2,z2],...[x3,y3,z3]]
#        list   unique(list)-> list                                                    # get non redundant list from a list
#        float  get_identity(seqA,seqB)                                                # get the identity pourcentage between seqA and seqB (Only the letters are saw). nb_ib between seqA and seqB / min(seqA,seqB)
#        list get_common_element(list1,list2)                                          # get a list with the common elements from 2 lists
#        float get_cover(master,seq)                                                   # get the coverage between a master sequence and a sequence (these sequence has to be align ! )
####################################################################

import math
import string


# Euclidian distance in 3D
def get_dist(atom1,atom2):
    """Give euclidean distance between 2 atoms"""
    dist=0.0
    for i in range(3):
        dist = dist+(atom1[i]-atom2[i])*(atom1[i]-atom2[i])
    return math.sqrt(dist)

def get_squared_dist(atom1,atom2):
    """Give euclidean distance between 2 atoms"""
    dist=0.0
    for i in range(3):
        dist = dist+(atom1[i]-atom2[i])*(atom1[i]-atom2[i])
    return dist


# Get a point from directionnal vector
def get_point(vector,ref_point,lambda_v):
    """Compute a point from a referent point, distance, and directionnal vector"""
    x = lambda_v * vector[0] + ref_point[0]
    y = lambda_v * vector[1] + ref_point[1]
    z = lambda_v * vector[2] + ref_point[2]
    return [x,y,z]



# Get colinear vector from 2 points
def get_vector(pa,pb):
    """Get colinear vector from pa to pb"""
    x=pb[0]-pa[0]
    y=pb[1]-pa[1]
    z=pb[2]-pa[2]
    lake = 1.0- get_dist(pa,pb)
    a= get_point([x,y,z],pb,lake / get_dist(pa,pb))

    xf=a[0]-pa[0]
    yf=a[1]-pa[1]
    zf=a[2]-pa[2]

    return [xf,yf,zf]



def get_min_dist(r1,r2):
    """ Get the minimal distance between 2 dictionnaries of atoms """
    min = 99999
    for i in r1.keys():
        for j in r2.keys():
            if i != "order" and j != "order":
                d = get_dist(r1[i],r2[j])
                if d < min :
                    min = d
    return min


def get_barycenter(coord):
    """ Give the barycenter of an ensemble of coordinates given in list coord=[[x1,y1,z1],[x2,y2,z2],...[x3,y3,z3]] """
    sumx = 0
    sumy = 0
    sumz = 0
    #print coord
    for i in coord.keys():
        sumx = sumx + coord[i][0]
        sumy = sumy + coord[i][1]
        sumz = sumz + coord[i][2]

    Xmean = float(sumx)/len(coord)
    Ymean = float(sumy)/len(coord)
    Zmean = float(sumz)/len(coord)
    bary = [Xmean,Ymean,Zmean]
    return bary

def get_barycenter_from_list(coord):
    """ Give the barycenter of an ensemble of coordinates given in list coord=[[x1,y1,z1],[x2,y2,z2],...[x3,y3,z3]] """
    sumx = 0
    sumy = 0
    sumz = 0
    #print coord
    for i in coord:
        sumx = sumx + i[0]
        sumy = sumy + i[1]
        sumz = sumz + i[2]

    Xmean = float(sumx)/len(coord)
    Ymean = float(sumy)/len(coord)
    Zmean = float(sumz)/len(coord)
    bary = [Xmean,Ymean,Zmean]
    return bary



def unique(s):
    """ Get unique element in a list"""
    se=[]
    for i in s:
        if not(i in se):
            se.append(i)
    return se

                                                
def get_identity(seqA,seqB):
    """ Compute the pourcentage of identity between two sequences """
    len_seqA = len(string.replace(seqA,"-",""))
    len_seqB = len(string.replace(seqB,"-",""))
    id=0
    if(len_seqA < len_seqB):
        min_len = len_seqA
    else:
        min_len = len_seqB
    for i in range(len(seqA)):
        if(seqA[i] == seqB[i] and seqA[i] != "-"):
            id = id + 1
    identity = float(float((id*100))/float(min_len))
    return identity
     

"""
def get_common_element(l1,l2):
    
    if min(len(l1),len(l2)) == l1:
        m1 = l1
        m2 = l2
    else:
        m1 = l2
        m2 = l1
    return [i for i in m1 if i in m2]
"""

def get_common_element(*lis):
    """
    Input : a list of list or several lists 
    Output : a list containing common elements between the multi lists
    Method uses the set() built-in is used so that intersection_update can be applied
    """ 
    if len(lis)==1:
        #
        # The input is in the form pmath.get_common_element([list1,list2,list3,...])
        #
        slis = [set(i) for i in lis[0]]
    else:
        #
        # The input is in the form pmath.get_common_element(list1,list2,list3,...)
        #
        slis = [set(i) for i in lis]
    [slis[0].intersection_update(i) for i in slis]
    return list(slis[0])

  

def get_ratio_length(master,seq):
    """ Get the cover from a master sequence """
    len_master = len(string.replace(master,"-",""))
    len_seq    = len(string.replace(seq,"-",""))

    return float(len_seq) * 100 / float(len_master)


def get_cover(seq_pdb,seq_cur):
    """    Given 2 aligned sequences (the pdb and an other one)      """
    """ compute the rate of aligned residues in the pdb sequence     """
    """ exemple:                                                     """
    """ pdb  ---PAVELETKARL----------AIVRY--MVFNARV--                """
    """ seq  GSTPAVILL-----TRYNRDEC---LVRF---------KS                """
    """ len_ali = 10, len_seq_pdb = 22, rat_ali_pdb = 10*100/22 = 45 """

    
    len_seq_pdb = len(string.replace(seq_cur,"-",""))
    len_seq_cur = len(string.replace(seq_cur,"-",""))
    len_ali = 0

    for i in range(len(seq_pdb)):
        if(seq_pdb[i] != "-" and seq_cur[i] != "-"):
            len_ali = len_ali + 1
    rat_ali_pdb = float(float((len_ali*100))/float(len_seq_pdb))

    return rat_ali_pdb
    """ End """



