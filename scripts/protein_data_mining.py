# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 09:47:34 2024

@author: Oscar Saiz Gutierrez
"""

# =====================================
# A
import requests  # Importo la librería 'requests'.

def descarga_pdb(pdb_ids): # Creo la función para descargar archivos PDBx/mmCIF a partir de una lista de IDs.
    for pdb_id in pdb_ids:
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"  # URL de descarga usando el ID de cada proteína.
        output_path = f"{pdb_id}.cif"  # Para crear el nombre del archivo donde guardaré cada archivo descargado.
        response = requests.get(url)  # Descargar el archivo.
        if response.status_code == 200:  # Verificar si la descarga fue bien.
            with open(output_path, "w") as file:  # Abrir el archivo en modo de escritura y guardar el contenido.
                file.write(response.text)
            print(f"File downloaded and saved: '{output_path}'.")  # Confirmar que la descarga fue exitosa.
        else:
            print(f"Error downloading file: {pdb_id}")  # Aviso si hay  error en la descarga.

# Lista de IDs de proteínas para descargar.
pdb_ids = ['1tup', '2xyz', '3def', '4ogq', '5jkl', '6mno', '7pqr', '8stu', '9vwx', '10yza']  # Las dos ultimas dan error, asi se comprueba que el aviso funciona

# Corro la función para descargar todos los archivos de la lista.
descarga_pdb(pdb_ids)

# =====================================
# B
import pandas as pd  # Importo la libreria pandas

# Creo la función para obtener el ID de UniProt a partir de un ID de PDB
def obtener_uniprot_id_desde_pdb(pdb_id):
    url = "https://rest.uniprot.org/idmapping/run"  # Dirección de la API de UniProt para convertir IDs
    datos = {'from': 'PDB', 'to': 'UniProtKB', 'ids': pdb_id}  
    respuesta = requests.post(url, data=datos)  
    if respuesta.status_code == 200:  
        job_id = respuesta.json()['jobId']  # Guardo el ID
        url_datos = f"https://rest.uniprot.org/idmapping/results/{job_id}"  # Dirección para los resultados
        resultado_final = requests.get(url_datos).json()  # Obtengo los resultados en formato JSON
        if 'results' in resultado_final and resultado_final['results']:  # Verifico si hay resultados
            uniprot_id = resultado_final['results'][0]['to']  # Obtengo el primer ID de UniProt
            print(f'ID de UniProt encontrado: {uniprot_id}')  # Muestro el ID de UniProt encontrado
            return uniprot_id  # ID de UniProt
        else:
            print(f"No hay  resultados para el ID de PDB {pdb_id}.")  # Mensaje si no hay resultados
            return None  # None si no se encuentra el ID de UniProt
    else:
        print(f"No se pudo obtener el ID de UniProt para el ID de PDB {pdb_id}.")  # Mensaje si falla
        return None  #  None en caso de error en la solicitud

# Función para obtener información de UniProt
def obtener_informacion_uniprot(uniprot_id):
    url = f"https://www.uniprot.org/uniprotkb/{uniprot_id}.json"  # Dirección de UniProt con el ID específico
    respuesta = requests.get(url)  
    if respuesta.status_code == 200:  
        datos = respuesta.json()  # Convierto la respuesta en un diccionario
        fecha_publicacion = datos['entryAudit']['firstPublicDate']  # Fecha de publicación 
        fecha_modificacion = datos['entryAudit']['lastAnnotationUpdateDate']  # Última modificación
        revisado = "Swiss-Prot" if datos['entryType'] else "Trembl"  # Revisado o no
        nombre_gen = datos['genes'][0]['geneName']['value']  # Nombre principal del gen
        sinonimos = ', '.join([syn['value'] for syn in datos['genes'][0].get('synonyms', [])])  # Sinónimos del gen
        organismo = datos['organism']['scientificName']  # Nombre científico del organismo
        nombre_proteina = datos['proteinDescription']['recommendedName']['fullName']['value']  # Nombre de la proteína
        secuencia_aa = datos['sequence']['value']  # Secuencia de aminoácidos
        pdb_ids = ', '.join([ref['id'] for ref in datos['uniProtKBCrossReferences'] if ref['database'] == 'PDB'])  # IDs de PDB relacionados
        return pd.Series(
            [uniprot_id, fecha_publicacion, fecha_modificacion, revisado, nombre_gen, sinonimos, organismo, nombre_proteina, secuencia_aa, pdb_ids],
            index=['Uniprot_id', 'Fecha_publicacion', 'Fecha_modificacion', 'Revisado', 'Nombre_del_gen', 'Sinónimos', 'Organismo', 'Nombre_proteina', 'Secuencia_aa', 'PDB_ids'])  
    else:
        print(f"No se pudo obtener la información de UniProt para el ID {uniprot_id}.")  # Mensaje  si hay error
        return pd.Series([uniprot_id] + [None] * 9)  # Serie en blanco si no se encuentra información

pdb_id = '1tup'  # Defino el ID para el que quiero obtener información
df = pd.DataFrame(columns=['Uniprot_id', 'Fecha_publicacion', 'Fecha_modificacion', 'Revisado', 'Nombre_del_gen', 'Sinónimos', 'Organismo', 'Nombre_proteina', 'Secuencia_aa', 'PDB_ids'])  # Creo un DataFrame vacío con las columnas

uniprot_id = obtener_uniprot_id_desde_pdb(pdb_id)  # Obtengo el ID de UniProt para el ID de PDB que he definido

if uniprot_id:  # Si el ID de UniProt es válido
    df.loc[0] = obtener_informacion_uniprot(uniprot_id)  # Agrego la información al DataFrame en la primera fila

print(df)  # Muestro el DataFrame con la información obtenida

# =====================================
# C
import pandas as pd  # Importo pandas

# Creo la función para obtener el nombre del cofactor de una entrada de UniProt usando su ID
def obtener_cofactor_uniprot(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()  # Obtengo los datos en formato JSON
        # Buscar la información sobre el cofactor
        for comment in data.get("comments", []):  
            if comment["commentType"] == "COFACTOR":  # Identifico si el comentario es de tipo "COFACTOR"
                return comment["cofactors"][0]["name"]  # Devuelvo el nombre del primer cofactor encontrado                
    return None  # Devuelvo None si no se encuentra el cofactor



# Creo la funcion para obtener información detallada del cofactor desde PubChem
def obtener_info_pubchem(cofactor_nombre):
    url_search = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cofactor_nombre}/cids/JSON"
    response = requests.get(url_search)
    if response.status_code == 200:
        cid = response.json()["IdentifierList"]["CID"][0]
        url_compound = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularWeight,InChI,InChIKey,IUPACName/JSON"
        response = requests.get(url_compound)  
        if response.status_code == 200:
            compound_data = response.json()["PropertyTable"]["Properties"][0]  # Obtengo los detalles del compuesto
            
            df_cofactor = pd.DataFrame([{ # Creo un DataFrame con la información obtenida sobre el cofactor
                "Compuesto": cofactor_nombre,  # Nombre del cofactor
                "Pubchem_id": cid,  # ID de PubChem
                "Peso_molecular": compound_data["MolecularWeight"],  # Peso molecular
                "Inchi": compound_data["InChI"],  # InChI
                "Inchikey": compound_data["InChIKey"],  # InChIKey
                "Iupac_name": compound_data["IUPACName"]  # Nombre IUPAC
            }])
            return df_cofactor  # DataFrame con los datos del cofactor
        else:
            print(f"No se pudo obtener información detallada para el CID {cid}.")  # Mensaje si falla

    return pd.DataFrame()  

cofactor_nombre = obtener_cofactor_uniprot("P04637")  # ID de UniProt para p53
print("Cofactor found:", cofactor_nombre)  # Print cofactor name found

if cofactor_nombre:  # Solo buscar en PubChem si se encontró un cofactor
    df_cofactor = obtener_info_pubchem(cofactor_nombre)  # Llamo a la función para obtener los datos y el DataFrame
    print(df_cofactor)  

# =====================================
# 2. Manipulación de datos biológicos
# A
from Bio import PDB  # Importo el módulo PDB de Biopython
import pandas as pd  # Importo pandas

parser = PDB.MMCIFParser(QUIET=True)  # Incio el parser
structure = parser.get_structure("protein", "4ogq.cif")  # Cargo la estructura de la proteína en la variable "structure".
heteroatoms = []  # Creo una lista vacía para almacenar los nombres de las heteromolecules encontradas (excluyendo el agua).
# Iterar cada modelo, cadena y residuo:
for model in structure:  
    for chain in model:  
        for residue in chain:
            # Exlcuir el agua ('HOH')
            if residue.id[0] != ' ' and residue.resname != 'HOH':  
                heteroatoms.append(residue.resname)  # Agrego el nombre del residuo a la lista de heteromolecules.


print("Heteromolecules encontradas (sin incluir agua):", heteroatoms)  # Muestra la lista de heteromolecules encontradas.

# =====================================
# B
from Bio.PDB.MMCIF2Dict import MMCIF2Dict  # Importo MMCIF2Dict 
import pandas as pd  # Importo pandas 

cif_dict = MMCIF2Dict("4ogq.cif")  # Lee el archivo .cif y lo convierte en un diccionario.

# Extraigo información sobre las heteromolecules:
nonpoly_names = cif_dict.get('_pdbx_entity_nonpoly.name', []) 
nonpoly_ids = cif_dict.get('_pdbx_entity_nonpoly.comp_id', [])

# Creo un DataFrame con la información:
df_heteroatoms = pd.DataFrame({  # Defino el DataFrame con dos columnas: 'Nombre' y 'ID_3Letras'.
    'Nombre': nonpoly_names,  
    'ID de tres letras': nonpoly_ids 
})

# Para excluir el agua
df_heteroatoms = df_heteroatoms[df_heteroatoms['ID de tres letras'] != 'HOH']

print("\nInformación de heteromolecules:")
print(df_heteroatoms)

# =====================================
# C
import requests
import pandas as pd  # Importo pandas
from rdkit import Chem  # Importo Chem de RDKit 
from rdkit.Chem import AllChem  # Importo AllChem de RDKit 

# Creo un DataFrame a partir de la lista de heteromolecules sin duplicados y con las columnas SMILES y SDF:
df_heteroatoms_sinduplicados = pd.DataFrame({'ID_3Letras': list(set(heteroatoms))})
df_heteroatoms_sinduplicados['SMILES'] = None
df_heteroatoms_sinduplicados['SDF'] = None

# Creo la funcion para obtener el SMILES desde PubChem
def obtener_smiles(nombre_heteromolecula):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{nombre_heteromolecula}/property/CanonicalSMILES/JSON"  # URL para obtener SMILES.
    response = requests.get(url)  
    if response.status_code == 200: 
        try:
            return response.json()['PropertyTable']['Properties'][0]['CanonicalSMILES']  
        except (KeyError, IndexError):
            print(f"No se encontró SMILES para {nombre_heteromolecula}.")  # Mensaje si no se encuentra SMILES.
            return None   
        
# Recorrer cada heteromolecule en el DataFrame para obtener SMILES y convertirlo en SDF
for index, row in df_heteroatoms_sinduplicados.iterrows():
    nombre = row['ID_3Letras']  # Obtengo el identificador de tres letras.
    smiles = obtener_smiles(nombre)  # Corro la función para obtener el SMILES.
    if smiles:  
        df_heteroatoms_sinduplicados.at[index, 'SMILES'] = smiles  # Almaceno el SMILES en el DataFrame.
        mol = Chem.MolFromSmiles(smiles)  # Convierto el SMILES a un objeto de molecule.
        if mol:  # Verifico si se creó la molecule.
            AllChem.Compute2DCoords(mol)  # Calculo las coordenadas 2D.
            sdf = Chem.MolToMolBlock(mol)  # Convierto la molecule a formato SDF.
            df_heteroatoms_sinduplicados.at[index, 'SDF'] = sdf  # Almaceno el SDF en el DataFrame.
        else:
            print(f"No se pudo convertir SMILES a molecule para {nombre}.")  # Mensaje si falla la conversión.

print(df_heteroatoms_sinduplicados)  # Imprimo el DataFrame con la información actualizada.



# =====================================
# D
df_heteroatoms_sdf = df_heteroatoms_sinduplicados.dropna(subset=['SDF']) # Para filtrar e incluir solo las moleculas que tienen SDF
# Creo archivo .sdf:
with Chem.SDWriter("heteromolecules.sdf") as sdf_writer:
    for index, row in df_heteroatoms_sinduplicados.iterrows():
        try: #Convertir el SMILES a un objeto molecula usando RDKit
            mol = Chem.MolFromSmiles(row['SMILES'])
            if mol:
                mol.SetProp("Molecular Weight", str(AllChem.CalcExactMolWt(mol)))
                sdf_writer.write(mol)
        except:
            pass

# Para validar que las moleculas son correctas y tienen el campo "Peso MOlecular"
supplier = Chem.SDMolSupplier('heteromolecules.sdf')
for mol in supplier:
    if mol is not None:
        if mol.HasProp("Molecular Weight"):
            print(Chem.MolToMolBlock(mol))
            print('-------------------')
            print("Molecular Weight:", mol.GetProp("Molecular Weight"))
        else:
            print("Error: 'Molecular_weight' no encontrado en la molecule.")
    else:
        print("Invalid molecule encontrada en el archivo SDF.")


