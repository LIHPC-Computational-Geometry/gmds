import subprocess

# Lire les noms de formes depuis le fichier shapes.txt
with open('listShapes.txt', 'r') as file:
    shapes = file.readlines()

# Supprimer les espaces et les sauts de ligne pour chaque forme
shapes = [shape.strip() for shape in shapes]

# Définir le nombre d'arguments que vous souhaitez passer à chaque exécutionqa
# Par exemple, passer 3 formes à chaque appel
batch_size = 1

# Lancer le programme C++ pour chaque groupe de formes
for i in range(0, len(shapes), batch_size):
    # Prenez un sous-ensemble des formes à envoyer à chaque exécution
    arg2 = shapes[i:i + batch_size]
    arg0="../../mctsblock/tst/data/params.json"
    arg1="/home/bourmaudp/Documents/PROJETS/devClass_Edges_Size_Path/blocking-learning-data/input/"+arg2+".vtk"



    # Lancer le programme C++ avec ces formes comme arguments
    print(f"Lancement du programme avec les formes: {', '.join(arg2)}")
    result = subprocess.run(['./main_blocker'] + arg0, capture_output=True, text=True)

    # Afficher la sortie du programme C++
    print(result.stdout)
    if result.stderr:
        print(f"Erreur : {result.stderr}")
