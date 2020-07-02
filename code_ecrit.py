import code_de_base as utils
import math
import numpy as np
import random
import scipy.stats as ss


############################################################ Projet3 ####################################################################
################################################## Analyse de séquences génomiques ######################################################
#########################################################################################################################################


# Merrouche Aymen
# Sidhoum Imad
# 2018/2019
# projet 3 3I005


#########################################################################################################################################
##########################################  Préliminaires : données et lecture des fichiers #############################################
#########################################################################################################################################


# nombre_nucleotide_total(nom_fichier)
# retourne le nombre total de nucleoitides dans un fichier qui represente un genome
# nom_fichier : le nom du fichicer fasta


def nombre_nucleotide_total(nom_fichier):
    # on recupere un dictionnaire qui represente chaque
    # sequence du fichier comme une lise d'entiers
    intermediate = utils.read_fasta(nom_fichier)
    cpt = 0
    # On récupère le nombre de nucléotides dans chaque séquence
    for sequence in intermediate:
        cpt += sum(utils.nucleotide_count(intermediate[sequence]))
    return cpt

# nombre_nucleotide_liste(nom_fichier)
# retourne le nombre total de nucleoitides selon le type de la nucleotides (liste à 4 éléments) dans un fichier qui represente un genome
# nom_fichier : le nom du fichicer fasta


def nombre_nucleotide_liste(nom_fichier):
    # on recupere un dictionnaire qui represente chaque
    # sequence du fichier comme une lise d'entiers
    intermediate = utils.read_fasta(nom_fichier)
    cpt = [0, 0, 0, 0]
    # On récupère le nombre de nucléotides dans chaque séquence
    for sequence in intermediate:
        cpt = np.array(cpt) + \
            np.array(utils.nucleotide_count(intermediate[sequence]))
    return cpt

# frequence_nucleotide(nom_fichier)
# retourne la frequence de chaque nucleotide dans un fichier qui represente un genome
# nom_fichier : le nom du fichicer fasta


"""
def frequence_nucleotide(nom_fichier):
    # on recupere un dictionnaire qui represente chaque
    # sequence du fichier comme une lise d'entiers
    intermediate = utils.read_fasta(nom_fichier)
    cpt = [0, 0, 0, 0]
    # On récupère les fréquences de nucléotides dans chaque séquence
    for sequence in intermediate:
        cpt += utils.nucleotide_frequency(intermediate[sequence])
    return cpt / len(intermediate)
"""


def frequence_nucleotide(nom_fichier):
    # on recupere un dictionnaire qui represente chaque
    # sequence du fichier comme une lise d'entiers
    intermediate = utils.read_fasta(nom_fichier)
    cpt = []
    # On récupère sous forme de liste les sequences concténées
    for sequence in intermediate:
        cpt += intermediate[sequence]
    # on renvoie les fréquences
    return utils.nucleotide_frequency(cpt)

# logproba(liste_entiers, m)
# calcule la log-probabilité d’une séquence étant donné les fréquences des lettres m
# liste_entiers : la séquence
# m : la fréquence des nucléotides


def logproba(liste_entiers, m):
    proba = 0
    # on applique la formule
    for nucl in liste_entiers:
        proba += math.log(m[nucl])
    return proba


# logproba(liste_entiers, m)
# calcule la log-probabilité d’une séquence étant donné les fréquences des lettres m, et le compte de nucleotides
# liste_comptage : le comptage des nucléotides
# m : la fréquence des nucléotides

def logproba_optimale(liste_comptage, m):
    proba = 0
    # = nb(nucleotides)*log(frequence(nucleotide))
    for i in range(len(liste_comptage)):
        proba += liste_comptage[i]*math.log(m[i])
    return proba


#########################################################################################################################################
#########################################################################################################################################


#########################################################################################################################################
################################################  Annotation des régions promoteurs #####################################################
#########################################################################################################################################


# concat_sequence(nom_fichier)
# rend sous forme de liste d'entiers les séquences concaténées de chaque fichier
# nom_fichier : nom du fichier fasta

def concat_sequence(nom_fichier):
    # on récupère sous forme de dictionnaire d'entiers les séaquences du fichier
    intermediate = utils.read_fasta(nom_fichier)
    cpt = []
    # on concténes toutes les séquences
    for sequence in intermediate:
        cpt += intermediate[sequence]
    return cpt

# code(m, k)
# renvoie pour un mot m de taille k son indice dans le tableau ordonné lexicographiquement
# m : le mot sous forme de liste d'entiers
# k : la longuer du mot


def code(m, k):
    taille = k-1
    indice = 0
    # l’écriture du mot en base 4
    for i in m:
        indice += i*(4 ** taille)
        taille -= 1
    return indice

# code_inverse(indice, k)
# connaissant un indice i et la longueur du mot k renvoie la séquence de longueur k correspondante
# indice : indice lexicographique du mot
# k : taille du mot


def code_inverse(indice, k):
    reste = indice
    taille = k-1
    liste = []
    div = 0
    # on récupère chque position dans l’écriture du mot en base 4
    while taille != -1:
        div = reste // (4**taille)
        reste = reste % (4**taille)
        liste.append(div)
        taille -= 1
    return liste
# compte_mot(sequence, k)
# compte le nombre d’occurrences pour tous les mots de taille k dans une séquence d’ADN
# Retourne un dictionnaire dont les cléfs sont les codes des mots et les valeurs sont les nombres d'occurences
# sequence : la séquence sur laquelle on va compter (liste d'entiers)
# k : la longuer des mots que nous allons compter


def compte_mot(sequence, k):
    dico = {}
    # on parcours la séquence
    for i in range(len(sequence) - k + 1):
        # on récupère le mot courant
        mot = sequence[i:i+k]
        # on récupère le code du mot courant
        code_mot = code(mot, k)
        # on met à jour le comptage
        if code_mot in dico:
            dico[code_mot] += 1
        else:
            dico[code_mot] = 1
    return dico

# comptage_attendu(frequneces, k, l)
# predit le nombre de mots attendus de taille k, dans une sequence de taille l, en se basant sur les frequences d'apparition des lettres
# Retourne un dictionnaire dont les cléfs sont les codes des mots et les valeurs sont les nombres d'occurences
# frequneces : les fréquences des nucléotides
# k : la longueur des mots que nous voulons compter
# l : la longuer de la séquence


def comptage_attendu(frequneces, k, l):
    # le nombre de mots de taille k dans une seqeunce de longuer l
    nombre_mot_attendu = l - k + 1
    dico = {}
    # on peroucrs tous les mots de taille k en utilisant leur codes lexicographiques
    for i in range(4**k):
        mot = code_inverse(i, k)
        # = proba du mot * nombre de mots de taille k dans une seqeunce de longuer l
        dico[i] = math.exp(logproba(mot, frequneces))*nombre_mot_attendu
    return dico

#########################################################################################################################################
#########################################################################################################################################


#########################################################################################################################################
################################################  Simulation de séquences aléatoires ####################################################
#########################################################################################################################################

# simule_sequence(lg, m)
# génère une séquence aléatoire de longueur lg d’une composition donnée m
# lg : longueur de la séquence
# m : proportion des nucléotides


def simule_sequence(lg, m):
    resultat = []
    # on genre autant de lettres que specifié dans m
    for i in range(len(m)):
        resultat += [i for j in range(int(lg*m[i]))]
    # on mélange le résultat
    random.shuffle(resultat)
    return resultat

# nb_occurence(sequence, mot)
# retourne le nombre d'occurence d'un mot dans une sequence
# sequence : la séquence (liste d'entiers)
# mot : le mot (liste d'entiers)


def nb_occurence(sequence, mot):
    k = len(mot)
    cpt = 0
    # our chaque mot de taille k dans la séquence
    for i in range(len(sequence) - k + 1):
        word = sequence[i:i+k]
        # si ce mot est égale à celui pour qui nous comptons les occurences on incrémente
        if (word == mot):
            cpt += 1
    return cpt

#########################################################################################################################################
#########################################################################################################################################


#########################################################################################################################################
###############################################  Modèles de dinucléotides et trinucléotides #############################################
#########################################################################################################################################

# simulation_M(sequence)
# Estime M à partir des comptages des mots de longueur 2.
# sequence : la sequence en liste d'entiers


def simulation_M(sequence):
    # on compte les mots de longueur 2
    comptage = compte_mot(sequence, 2)
    c = 0
    # on initialise notre matrice
    M = [[0, 0, 0, 0] for i in range(4)]
    N = [0, 0, 0, 0]
    # pour chaque dinucléotide
    # on recupére le nombre total d'apparitions de cette dinucliotide
    # et on divise par le nombre de dicnucliotide qui commence par i car P(x/y) est une proba
    for i in [0, 1, 2, 3]:
        for j in [0, 1, 2, 3]:
            mot = [i, j]
            c = code(mot, 2)
            if c in comptage:
                M[i][j] = comptage[c]
                N[i] += comptage[c]
            else:
                M[i][j] = 0
    for i in [0, 1, 2, 3]:
        M[i] = np.array(M[i]) / N[i]
    return M

# simulation_dinucleotide(M, frequence, l)
# simule une séquence de longueur l avec le modèle de dinucléotides
# M : matrice de transitions
# frequence : fréquence des nucleotides
# l : longueur de la séquence en sortie


def simulation_dinucleotide(M, frequence, l):
    liste = []
    current_frequence = frequence
    for i in range(l):
        # on génére une lettre suivant une certaien proportion de nucléotides
        # cette proportion est la ligne dans M qui correspond à la lettre courante
        choix = np.random.choice([0, 1, 2, 3], 1, p=current_frequence)[0]
        liste.append(choix)
        current_frequence = M[choix]
    return liste
# proba_markov(m, frequences, matrice)
# Calcule la proba d'un mot dans une sequence
# m : le mot
# frequences : les fréquences de nucléotides
# matrice : la matrice de transitions


def proba_markov(m, frequences, matrice):
    prod = frequences[m[0]]
    # Il suffit d'utiliser la décomposition de la loi jointe dans une chaine de Markov.
    for i in range(1, len(m)):
        prod = prod*matrice[m[i-1]][m[i]]
    return prod

# comptage_attendu_markov(frequneces, matrice, k, l)
# predit le nombre de mots attendus de taille k, dans une sequence de taille l, en se basant sur les frequences d'apparition des lettres
# et la matrice de transitions
# Retourne un dictionnaire dont les cléfs sont les codes des mots et les valeurs sont les nombres d'occurences
# frequneces : les fréquences des nucléotides
# matrice : la matrice de transition
# k : la longueur des mots que nous voulons compter
# l : la longuer de la séquence


def comptage_attendu_markov(frequneces, matrice, k, l):
    # = proba du mot * nombre de mots de taille k
    nombre_mot_attendu = l - k + 1
    dico = {}
    # on peroucrs tous les mots de taille k
    for i in range(4**k):
        mot = code_inverse(i, k)
        # = proba du mot * nombre de mots de taille k
        dico[i] = proba_markov(mot, frequneces, matrice)*nombre_mot_attendu
    return dico


# detect_mot_chevauchant(m)
# renvoie si oui ou non le mot est chevauchant avec lui meme
# m : le mot
def detect_mot_chevauchant(m):
    dico = compte_mot(m, 2)
    if m[0] == m[len(m)-1]:
        return True
    for mot in dico:
        # si il y a un dinucleotide qui se repete
        if dico[mot] >= 2:
            return True
    return False


#########################################################################################################################################
#########################################################################################################################################


#########################################################################################################################################
##########################################################  Probabilités de mots ########################################################
#########################################################################################################################################

# estime_proba_empirique(dico, l, k, n)
# Fonction qui, à partir de la liste des comptages des mots de taille k et de la longueur de la séquence, calcul leur probabilité empirique
# dico : comptage des mots de taille k
# l : longueur des la chaine
# k : longueur des mots
# n : paramétre de la probabilité empirique P(x>=n)

def estime_proba_empirique(dico, l, k, n):
    proba = {}
    p = 0
    # on itére sur tout les mots de longuer k
    for word in range(4**k):
        proba[word] = 0
        for i in range(n):
            # on récupère le comptage du mot
            if word in dico:
                p = dico[word]
            else:
                p = 0
            # on calcule la proba
            proba[word] += ss.poisson.pmf(i, (l-k+1)*(p/(l-k+1)), loc=0)
        proba[word] = 1-proba[word]
    return proba

# estime_proba_empirique(dico, l, k, n)
# Fonction qui, à partir de la liste des comptages des mots de taille k et de la longueur de la séquence, calcul leur probabilité empirique
# l : longueur des la chaine
# k : longueur des mots
# n : paramétre de la probabilité empirique P(x>=n)
# frequences : frequences des nucleotides


def estime_proba_empirique_naive(l, k, n,frequences):
    proba = {}
    p = 0
    # on itére sur tout les mots de longuer k
    for word in range(4**k):
        proba[word] = 0
        # on calcule la proba
        probabi = math.exp(logproba(code_inverse(word,k), frequences))
        for i in range(n):
            proba[word] += ss.poisson.pmf(i, (l-k+1)*probabi, loc=0)
        proba[word] = 1-proba[word]
    return proba


#########################################################################################################################################
#########################################################################################################################################
