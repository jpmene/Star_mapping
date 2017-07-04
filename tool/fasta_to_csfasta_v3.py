#!/usr/bin/python
"""
l'outil permet la convertion d'un fichier fasta en un fichier csfasta avec le codage
du format color .

exemple : 

ACCTGNTAT
|||||| || 
A1021.133
"""


import sys
import time 
"""
fonction permettant de verifier quel couple de nucleotide correspond au code color 

input : 
	- i correspond a une position de la sequence 
	- j correspond a une autre position de la sequence (ex j = i+1)
	- sequence correspond a la sequence de nucleotide 
"""
def Condition(i , j,sequence ):
	read = ""
	if sequence[i].upper() == "A" :
		if sequence[j].upper() == "A":
			read = "0"
		elif sequence[j].upper() == "C":
			read = "1"
		elif sequence[j].upper() == "G":
			read = "2"
		elif sequence[j].upper() == "T":		
			read = "3"
		elif sequence[j].upper() == "N":		
			read = "."

	elif sequence[i].upper() =="C" :
		if sequence[j].upper() == "A":
			read = "1"
		elif sequence[j].upper() == "C":
			read = "0"
		elif sequence[j].upper() == "G":
			read = "3"
		elif sequence[j].upper() == "T":		
			read = "2"
		elif sequence[j].upper() == "N":		
			read = "."

	elif sequence[i].upper() =="G" :
		if sequence[j].upper() == "A":
			read = "2"
		elif sequence[j].upper() == "C":
			read = "3"
		elif sequence[j].upper() == "G":
			read = "0"
		elif sequence[j].upper() == "T":		
			read = "1"
		elif sequence[j].upper() == "N":		
			read = "."

	elif sequence[i].upper() =="T": 
		if sequence[j].upper() == "A":
			read = "3"
		elif sequence[j].upper() == "C":
			read = "2"
		elif sequence[j].upper() == "G":
			read = "1"
		elif sequence[j].upper() == "T":		
			read = "0"
		elif sequence[j].upper() == "N":		
			read = "."

	elif sequence[i].upper() =="N": 
		if sequence[j].upper() == "A":
			read = "."
		elif sequence[j].upper() == "C":
			read = "."
		elif sequence[j].upper() == "G":
			read = "."
		elif sequence[j].upper() == "T":		
			read = "."
		elif sequence[j].upper() == "N":		
			read = "."


	return read


fichier = open((sys.argv[1]), "r")
csfasta = open((sys.argv[2]), "w")

# lecture du fichier par lignes 
line = fichier.readlines()


tps1 = time.clock() 

# pour chaque ligne 
for sequence in line:
	# recupere la premier ligne avec ">"
	if sequence[0] == ">":
		csfasta.write(sequence) 
		end =""
		a = True 
	else:
		# enregistre le premier nucleotide de la sequence 
		if a == True :
			csfasta.write(sequence[0].upper())
			a= False 

		read =''

		
		# ajout le dernier nucleotide de la sequence precedente 
		sequence = end+sequence
		#print(sequence)
		# converti les nucleotide en color 
		for i in range(0,len(sequence)-1):
			#print(read)
			read += Condition(i,i+1,sequence)


		# recupere le dernier nucleotide de la sequence 			
		end = sequence[-2]
		
		#print(read)		
		# ecrit la ligne dans le nouveaux fichier 
		csfasta.write( read+"\n")

fichier.close()
csfasta.close()

tps2 = time.clock() 
print(tps2 - tps1)
